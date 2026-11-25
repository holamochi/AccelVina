#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import os
import re
import resource
import shutil
import socket
import statistics
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


KERNELS = ("Kernel1_Opt.bin", "Kernel2_Opt.bin")
SCORE_ROW = re.compile(r"\s*\d+\s+-?\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s*")


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def case_args(case, suite, input_root, batch_root):
    arguments = ["--receptor", str(input_root / case["receptor"])]
    for axis, center, size in zip("xyz", case["center"], case["size"], strict=True):
        arguments += [f"--center_{axis}", str(center), f"--size_{axis}", str(size)]
    settings = suite["settings"] | case.get("settings", {})
    for option, value in settings.items():
        arguments += [f"--{option}", str(value)]
    if "flex" in case:
        arguments += ["--flex", str(input_root / case["flex"])]
    if "ligand_directory" not in case:
        arguments += ["--ligand", str(input_root / case["ligand"]), "--out", "out/docked.pdbqt"]
        return arguments, {"docked.pdbqt"}
    directory = input_root / case["ligand_directory"]
    ordered = batch_root / case["id"]
    ordered.mkdir()
    # Upstream consumes directory order and uses one random stream for the entire batch.
    for name in reversed(case["ligand_order"]):
        shutil.copyfile(directory / name, ordered / name)
    if [p.name for p in ordered.iterdir()] != case["ligand_order"]:
        raise ValueError(f"{case['id']}: tmpfs did not preserve the frozen batch order")
    arguments += ["--ligand_directory", str(ordered), "--output_directory", "out"]
    return arguments, {f"{Path(name).stem}_out.pdbqt" for name in case["ligand_order"]}


def read_outputs(out_dir, expected, stdout, stderr, num_modes):
    files = {path.name: path for path in out_dir.iterdir()}
    if set(files) != expected:
        raise ValueError(f"PDBQT output names differ (missing {sorted(expected - files.keys())}, extra {sorted(files.keys() - expected)})")

    outputs = {}
    for name, path in sorted(files.items()):
        lines = path.read_text().splitlines()
        starts = [i for i, line in enumerate(lines) if line.startswith("MODEL ")]
        ends = [i for i, line in enumerate(lines) if line.startswith("ENDMDL")]
        if not 1 <= len(starts) == len(ends) <= num_modes:
            raise ValueError(f"{name}: missing or malformed models")
        scores = []
        for start, end in zip(starts, ends, strict=True):
            model = lines[start:end]
            remarks = [line.split()[3:] for line in model if line.startswith("REMARK VINA RESULT:")]
            if len(remarks) != 1 or len(remarks[0]) != 3 or not all(math.isfinite(float(v)) for v in remarks[0]):
                raise ValueError(f"{name}: invalid Vina result row")
            if not any(line.startswith(("ATOM  ", "HETATM")) for line in model):
                raise ValueError(f"{name}: empty model")
            scores.append(float(remarks[0][0]))
        outputs[name] = {"sha256": sha256(path), "modes": len(scores), "affinity_kcal_mol": scores}

    score_rows = [" ".join(line.split()) for line in stdout.read_text().splitlines() if SCORE_ROW.fullmatch(line)]
    if len(score_rows) != sum(output["modes"] for output in outputs.values()):
        raise ValueError("stdout pose table is incomplete")
    if stderr.stat().st_size:
        raise ValueError("stderr is nonempty")
    return {"pdbqt": outputs, "stdout_pose_rows": score_rows}


def run_case(command, expected, work, env, timeout, num_modes):
    (work / "out").mkdir(parents=True)
    stdout, stderr = work / "stdout.log", work / "stderr.log"
    record = {"argv": command, "directory": str(work)}
    with stdout.open("wb") as out_log, stderr.open("wb") as err_log:
        start = time.perf_counter()
        try:
            result = subprocess.run(command, cwd=work, env=env, stdout=out_log, stderr=err_log, timeout=timeout, check=False)
            record["exit_code"] = result.returncode
        except subprocess.TimeoutExpired:
            record.update(exit_code=None, error=f"Timed out after {timeout} seconds")
        record["wall_s"] = time.perf_counter() - start
    if record["exit_code"] == 0:
        try:
            record["signature"] = read_outputs(work / "out", expected, stdout, stderr, num_modes)
        except (OSError, ValueError) as exc:
            record["error"] = str(exc)
    elif record["exit_code"] is not None:
        record["error"] = f"Exit code {record['exit_code']}"
    return record


def main():
    parser = argparse.ArgumentParser(description="Run a frozen AutoDock-Vina-GPU 2.1 suite on one NVIDIA GPU.")
    parser.add_argument("--suite", type=Path, default=Path(__file__).parent / "demo_test" / "demo_test.json")
    parser.add_argument("--baseline-runtime", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New run directory")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--candidate-kernels", type=Path)
    parser.add_argument("--reference", type=Path, help="Frozen base_result JSON for this suite")
    args = parser.parse_args()
    if args.repeats < 1 or (args.candidate_kernels and not args.candidate) or (args.candidate and not args.reference):
        parser.error("Use positive repeats, and supply --reference when running a candidate.")
    with tempfile.TemporaryDirectory(dir="/dev/shm", prefix="vina-gpu-batch-") as tmp:
        return benchmark(args, Path(tmp))


def benchmark(args, batch_root):
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    suite_path = args.suite.resolve(strict=True)
    suite = json.loads(suite_path.read_text())
    runtime, output = args.baseline_runtime.resolve(strict=True), args.output.resolve()
    input_root = suite_path.parent
    for name, expected in suite["inputs_sha256"].items():
        if sha256(input_root / name) != expected:
            raise ValueError(f"Input SHA-256 differs: {name}")
    fingerprint = {}
    for line in (runtime / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if sha256(runtime / name) != expected:
            raise ValueError(f"Baseline runtime SHA-256 differs: {name}")
        fingerprint[name] = expected

    gpu = subprocess.check_output(["nvidia-smi", "--query-gpu=uuid,name,driver_version", "--format=csv,noheader"], text=True).splitlines()
    if len(gpu) != 1:
        raise ValueError("Expose exactly one NVIDIA GPU.")
    busy = subprocess.check_output(["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], text=True).strip()
    if busy:
        raise ValueError("The GPU already has active compute processes.")
    cpus = suite["cpu_affinity"]
    os.sched_setaffinity(0, cpus)
    env = {
        "PATH": os.environ["PATH"],
        "LC_ALL": "C",
        "LD_LIBRARY_PATH": f"{runtime / 'lib'}:/usr/lib/x86_64-linux-gnu",
        "OCL_ICD_VENDORS": str(runtime / "icd"),
        "CUDA_CACHE_DISABLE": "0",
    }
    report = {
        "suite": suite["name"], "suite_sha256": sha256(suite_path), "baseline_sha256": fingerprint,
        "hostname": socket.gethostname(), "gpu": gpu[0], "cpu_affinity": cpus,
        "started_utc": datetime.now(timezone.utc).isoformat(), "repeats": args.repeats,
        "runtime": str(runtime), "cases": {},
    }
    reference = json.loads(args.reference.read_text()) if args.reference else None
    if reference and (not reference["valid"] or any(reference[key] != report[key]
                      for key in ("suite_sha256", "baseline_sha256", "gpu", "cpu_affinity"))):
        raise ValueError("Reference does not match the suite, runtime, GPU and CPU affinity.")
    versions = [("baseline", runtime / "AutoDock-Vina-GPU-2-1", runtime)]
    if args.candidate:
        candidate = args.candidate.resolve(strict=True)
        kernels = args.candidate_kernels.resolve(strict=True) if args.candidate_kernels else runtime
        versions.append(("candidate", candidate, kernels))
        report["candidate_sha256"] = {"entrypoint": sha256(candidate)} | {name: sha256(kernels / name) for name in KERNELS}

    cases = [(case, *case_args(case, suite, input_root, batch_root)) for case in suite["cases"]]
    output.mkdir()
    cache_paths = {name: output / "driver-cache" / name for name, _, _ in versions}
    for cache in cache_paths.values():
        cache.mkdir(parents=True)
    for case, _, _ in cases:
        report["cases"][case["id"]] = {"runs": {name: [] for name, _, _ in versions}}

    for repeat in range(1, args.repeats + 1):
        order = versions if repeat % 2 else versions[::-1]
        for case, cli, expected in cases:
            group = report["cases"][case["id"]]
            for name, executable, kernels in order:
                work = output / f"repeat-{repeat}" / case["id"] / name
                command = [str(executable), *cli, "--opencl_binary_path", str(kernels)]
                num_modes = case.get("settings", {}).get("num_modes", suite["settings"]["num_modes"])
                run_env = env | {"CUDA_CACHE_PATH": str(cache_paths[name]), "TMPDIR": str(work)}
                record = run_case(command, expected, work, run_env, suite["timeout_s"], num_modes)
                group["runs"][name].append(record)
                print(f"{repeat}/{args.repeats} {case['id']} {name} {record['wall_s']:.3f}s {record.get('error', 'completed')}", flush=True)

    for case_id, group in report["cases"].items():
        first = group["runs"]["baseline"][0]
        expected = reference["cases"][case_id]["reference"] if reference else first.get("signature")
        group["reference"] = expected
        for runs in group["runs"].values():
            for run in runs:
                run["valid"] = expected is not None and run.get("signature") == expected
                if not run["valid"]:
                    run.setdefault("error", "PDBQT or stdout pose table differs from the reference")
        group["valid"] = all(run["valid"] for runs in group["runs"].values() for run in runs)
        group["median_wall_s"] = {
            name: statistics.median(run["wall_s"] for run in runs) for name, runs in group["runs"].items()
        } if group["valid"] else {}

    report["valid"] = all(group["valid"] for group in report["cases"].values())
    if report["valid"]:
        report["total_median_wall_s"] = {
            name: sum(group["median_wall_s"][name] for group in report["cases"].values()) for name, _, _ in versions
        }
        if args.candidate:
            totals = report["total_median_wall_s"]
            report["speedup"] = totals["baseline"] / totals["candidate"]
    (output / "result.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"{output / 'result.json'} (valid={report['valid']})", flush=True)
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
