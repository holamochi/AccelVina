#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import os
import re
import resource
import shlex
import signal
import socket
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


GPU_EXAMPLE = Path(__file__).resolve().parents[2]
sys.dont_write_bytecode = True
from isolated import isolated_command  # noqa: E402

sys.path.insert(0, str(GPU_EXAMPLE))
from run_benchmark import KERNELS, SCORE_ROW, case_args  # noqa: E402


DEFAULT_SUITE = GPU_EXAMPLE / "demo_test" / "demo_test.json"
DEVICE_ENV = ("MUSA_VISIBLE_DEVICES", "CUDA_VISIBLE_DEVICES", "OCL_ICD_VENDORS",
              "OCL_ICD_FILENAMES", "OPENCL_VENDOR_PATH", "GPU_DEVICE_ORDINAL", "VINA_EXPECT_PCI")
STREAM_ENV = ("PATH", "LD_LIBRARY_PATH", "LD_PRELOAD")
DEVICE_ROW = re.compile(r"^GPU Device:[ \t]*(.+)$", re.MULTILINE)
RESERVED_OPTIONS = {"receptor", "ligand", "ligand_directory", "flex", "out", "output_directory",
                    "opencl_binary_path", "log", "config"} | {
    f"{box}_{axis}" for box in ("center", "size") for axis in "xyz"
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_hashes(root, fingerprints, kind):
    for name, expected in fingerprints.items():
        path = root / name
        if Path(name).is_absolute() or ".." in Path(name).parts or not path.resolve().is_relative_to(root):
            raise ValueError(f"{kind} escapes its directory: {name}")
        if sha256(path) != expected:
            raise ValueError(f"{kind} SHA-256 differs: {name}")


def load_suite(suite_path, batch_root):
    suite = json.loads(suite_path.read_text())
    if suite["schema_version"] != 1 or not suite["cases"] or suite["timeout_s"] <= 0:
        raise ValueError("Expected a nonempty schema_version 1 suite with a positive timeout")
    verify_hashes(suite_path.parent, suite["inputs_sha256"], "Input")
    cases = {}
    for case in suite["cases"]:
        case_id = case["id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", case_id) or case_id in cases:
            raise ValueError(f"Invalid or duplicate case ID: {case_id}")
        inputs = [case[key] for key in ("receptor", "flex") if key in case]
        if ("ligand" in case) == ("ligand_directory" in case):
            raise ValueError(f"{case_id}: expected exactly one ligand or ligand_directory")
        if "ligand_directory" in case:
            order = case["ligand_order"]
            if not order or len(set(order)) != len(order) or any(
                Path(name).name != name or not name.endswith(".pdbqt") for name in order
            ):
                raise ValueError(f"{case_id}: invalid ligand_order")
            inputs.extend(str(Path(case["ligand_directory"]) / name) for name in order)
            outputs = [f"{Path(name).stem}_out.pdbqt" for name in order]
        else:
            inputs.append(case["ligand"])
            outputs = ["docked.pdbqt"]
        unhashed = set(inputs) - suite["inputs_sha256"].keys()
        if unhashed:
            raise ValueError(f"{case_id}: inputs missing SHA-256 values: {sorted(unhashed)}")
        settings = suite["settings"] | case.get("settings", {})
        if RESERVED_OPTIONS & settings.keys():
            raise ValueError(f"{case_id}: settings override the benchmark input or output paths")
        if type(settings["num_modes"]) is not int or settings["num_modes"] < 1:
            raise ValueError(f"{case_id}: num_modes must be a positive integer")
        command, _ = case_args(case, suite, suite_path.parent, batch_root)
        cases[case_id] = (command, outputs, settings["num_modes"])
    return suite, cases


def read_signature(out_dir, expected, stdout, stderr, num_modes):
    if stderr.stat().st_size:
        raise ValueError("stderr is nonempty")
    actual = {path.name for path in out_dir.iterdir()}
    if actual != set(expected):
        raise ValueError(f"PDBQT names differ (missing {sorted(set(expected) - actual)}, "
                         f"extra {sorted(actual - set(expected))})")
    outputs, file_rows = {}, []
    for name in expected:
        path = out_dir / name
        lines = path.read_text().splitlines()
        starts = [i for i, line in enumerate(lines) if line.startswith("MODEL ")]
        ends = [i for i, line in enumerate(lines) if line.startswith("ENDMDL")]
        if not 1 <= len(starts) == len(ends) <= num_modes:
            raise ValueError(f"{name}: missing or malformed models")
        affinities = []
        for mode, (start, end) in enumerate(zip(starts, ends, strict=True), 1):
            if start >= end or (mode > 1 and start <= ends[mode - 2]) or lines[start].split() != ["MODEL", str(mode)]:
                raise ValueError(f"{name}: invalid MODEL order")
            model = lines[start + 1:end]
            remarks = [line.split()[3:] for line in model if line.startswith("REMARK VINA RESULT:")]
            if len(remarks) != 1 or len(remarks[0]) != 3:
                raise ValueError(f"{name}: missing Vina result row in model {mode}")
            values = tuple(map(float, remarks[0]))
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{name}: non-finite Vina result row in model {mode}")
            if not any(line.startswith(("ATOM  ", "HETATM")) for line in model):
                raise ValueError(f"{name}: empty model {mode}")
            file_rows.append((mode, *values))
            affinities.append(values[0])
        outputs[name] = {"sha256": sha256(path), "modes": len(affinities), "affinity_kcal_mol": affinities}
    score_rows = [" ".join(line.split()) for line in stdout.splitlines() if SCORE_ROW.fullmatch(line)]
    stdout_rows = []
    for row in score_rows:
        mode, *values = row.split()
        stdout_rows.append((int(mode), *map(float, values)))
    if stdout_rows != file_rows:
        raise ValueError("stdout pose table differs from PDBQT pose rows in frozen ligand order")
    return {"pdbqt": outputs, "stdout_pose_rows": score_rows}


def run_case(command, expected, work, env, timeout, num_modes, expected_device_name, isolation=None):
    (work / "out").mkdir(parents=True)
    stdout, stderr = work / "stdout.log", work / "stderr.log"
    record = {"argv": command, "directory": str(work), "exit_code": None}
    with stdout.open("wb") as out_log, stderr.open("wb") as err_log:
        start = None
        try:
            argv = command
            process_env = env | {"TMPDIR": str(work)}
            stdin = None
            if isolation:
                library_dirs = list(isolation["library_dirs"])
                executable_lib = Path(command[0]).parent / "lib"
                if executable_lib.is_dir() and executable_lib not in library_dirs:
                    library_dirs.append(executable_lib)
                argv = isolated_command(command, work, env, **(isolation | {"library_dirs": library_dirs}))
                record["sandbox_argv"] = argv
                process_env = {}
                stdin = subprocess.DEVNULL
            start = time.perf_counter()
            with subprocess.Popen(argv, cwd=work, env=process_env, stdin=stdin, stdout=out_log,
                                  stderr=err_log, start_new_session=True) as process:
                try:
                    record["exit_code"] = process.wait(timeout=timeout)
                finally:
                    if process.returncode is None:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        process.wait()
        except subprocess.TimeoutExpired:
            record["error"] = f"Timed out after {timeout} seconds"
        except (OSError, ValueError) as exc:
            record["error"] = str(exc)
        record["wall_s"] = time.perf_counter() - start if start is not None else 0.0
    if record["exit_code"] == 0:
        try:
            log = stdout.read_text()
            devices = [name.strip() for name in DEVICE_ROW.findall(log)]
            record["reported_devices"] = devices
            if len(devices) > 1 or (expected_device_name and devices != [expected_device_name]):
                raise ValueError(f"Unexpected GPU Device records: {devices}")
            record["signature"] = read_signature(work / "out", expected, log, stderr, num_modes)
        except (OSError, ValueError) as exc:
            record["error"] = str(exc)
    elif record["exit_code"] is not None:
        record["error"] = f"Exit code {record['exit_code']}"
        if isolation:
            with stderr.open("rb") as error_log:
                message = error_log.read(2048).decode(errors="replace").strip()
            if message:
                record["error"] += f" ({message.splitlines()[0]})"
    return record


def version_artifacts(executable, backend, kernels):
    executable = executable.resolve(strict=True)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"Not an executable file: {executable}")
    paths = {"entrypoint": executable}
    if backend == "opencl":
        kernels = kernels.resolve(strict=True) if kernels else executable.parent
        paths.update({name: kernels / name for name in KERNELS})
    elif kernels:
        raise ValueError("OpenCL kernel directories apply only to versions with the opencl backend")
    version = {
        "backend": backend,
        "artifacts": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
        "prefix": [str(executable)],
        "suffix": ["--opencl_binary_path", str(kernels)] if kernels else [],
    }
    manifest = executable.parent / "SHA256SUMS"
    if manifest.is_file():
        fingerprints = {}
        for line in manifest.read_text().splitlines():
            checksum, name = line.split(maxsplit=1)
            fingerprints[name.removeprefix("*")] = checksum
        verify_hashes(manifest.parent, fingerprints, "Build artifact")
        version["manifest"] = {"path": str(manifest), "sha256": sha256(manifest), "files": fingerprints}
    return version


def version_environment(env, cache, extra):
    paths = {"HOME": cache / "home", "XDG_CACHE_HOME": cache / "xdg"}
    if extra:
        paths.update({name: cache / "native" for name in extra})
    for path in set(paths.values()):
        path.mkdir(parents=True)
    return env | {name: str(path) for name, path in paths.items()}


def probe_device(command, env):
    argv = shlex.split(command)
    result = subprocess.run(argv, env=env, capture_output=True, text=True, timeout=30, check=True)
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"Device probe must report exactly one visible GPU, got {len(lines)} lines")
    return {"argv": argv, "visible_gpu": lines[0], "stderr": result.stderr}


def load_a100_reference(path, suite_digest, cases):
    path = path.resolve(strict=True)
    reference = json.loads(path.read_text())
    if (not reference["valid"] or not re.search(r"\bA100\b", reference["gpu"]) or
            reference["suite_sha256"] != suite_digest or reference["cases"].keys() != cases.keys()):
        raise ValueError("A100 reference must be valid and match the suite SHA-256 and case IDs")
    signatures = {case_id: group["reference"] for case_id, group in reference["cases"].items()}
    if any(not group["valid"] or not signatures[case_id] for case_id, group in reference["cases"].items()):
        raise ValueError("A100 reference is missing a valid signature")
    return {"path": str(path), "sha256": sha256(path), "gpu": reference["gpu"]}, signatures


def compare_signatures(actual, expected):
    if actual is None:
        return {"status": "unavailable"}
    files = actual["pdbqt"].keys() | expected["pdbqt"].keys()
    differing = sorted(name for name in files if actual["pdbqt"].get(name) != expected["pdbqt"].get(name))
    stdout_matches = actual["stdout_pose_rows"] == expected["stdout_pose_rows"]
    return {"status": "match" if stdout_matches and not differing else "different",
            "pdbqt_differences": differing, "stdout_pose_rows_match": stdout_matches}


def summarize(report, versions):
    for group in report["cases"].values():
        group["valid"] = all(run["valid"] for runs in group["runs"].values() for run in runs)
        group["median_wall_s"] = {
            version: statistics.median(run["wall_s"] for run in runs) for version, runs in group["runs"].items()
        } if group["valid"] else {}
    if "a100_reference" in report:
        report["a100_reference"]["run_status_counts"] = {
            name: dict(Counter(run["a100_comparison"]["status"] for group in report["cases"].values()
                               for run in group["runs"][name])) for name in versions
        }
    report["valid"] = all(group["valid"] for group in report["cases"].values())
    if report["valid"]:
        totals = {name: sum(group["median_wall_s"][name] for group in report["cases"].values()) for name in versions}
        report["total_median_wall_s"] = totals
        report["speedup"] = totals["baseline"] / totals["candidate"]


def benchmark(args, batch_root):
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    suite_path = args.suite.resolve(strict=True)
    suite_digest = sha256(suite_path)
    suite, cases = load_suite(suite_path, batch_root)
    if args.cpu_affinity:
        requested = {int(cpu) for cpu in args.cpu_affinity.split(",")}
        available = os.sched_getaffinity(0)
        if not requested <= available:
            raise ValueError(f"Requested CPU IDs outside the current affinity: {sorted(requested - available)}")
        os.sched_setaffinity(0, requested)
    versions = {
        "baseline": version_artifacts(args.baseline, args.backend, args.baseline_kernels),
        "candidate": version_artifacts(args.candidate, args.candidate_backend or args.backend, args.candidate_kernels),
    }
    env = os.environ.copy() | {"LC_ALL": "C"}
    report = {
        "schema_version": 1,
        "suite": suite["name"], "suite_sha256": suite_digest, "gpu": args.gpu,
        "hostname": socket.gethostname(), "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "started_utc": datetime.now(timezone.utc).isoformat(), "repeats": args.repeats,
        "runtime_env": {name: env[name] for name in (*DEVICE_ENV, *STREAM_ENV, "VINA_MUSA_PROFILE") if name in env},
        "versions": versions,
        "cases": {case_id: {"runs": {name: [] for name in versions}} for case_id in cases},
    }
    if args.device_probe:
        report["device_probe"] = probe_device(args.device_probe, env)
    if args.expected_device_name:
        report["expected_device_name"] = args.expected_device_name
    isolation = None
    if args.isolate:
        isolation = {"bwrap": str(args.bwrap or "/usr/bin/bwrap"), "device_nodes": args.device_node,
                     "library_dirs": args.library_dir, "cache_env": args.cache_env}
        report["isolation"] = isolation
    a100_signatures = None
    if args.reference:
        report["a100_reference"], a100_signatures = load_a100_reference(args.reference, suite_digest, cases)

    output = args.output.resolve()
    output.mkdir()
    run_env = {name: version_environment(env, output / "driver-cache" / name, args.cache_env) for name in versions}
    for name, settings in versions.items():
        settings["cache_env"] = {key: run_env[name][key] for key in ("HOME", "XDG_CACHE_HOME", *args.cache_env)}
    report["repeat_order"] = []
    device_labels = {}
    for repeat in range(1, args.repeats + 1):
        order = list(versions) if repeat % 2 else list(reversed(versions))
        report["repeat_order"].append(order)
        for case_id, (cli, expected, num_modes) in cases.items():
            group = report["cases"][case_id]
            for name in order:
                settings = versions[name]
                work = output / f"repeat-{repeat}" / case_id / name
                command = [*settings["prefix"], *cli, *settings["suffix"]]
                run = run_case(command, expected, work, run_env[name], suite["timeout_s"], num_modes,
                               args.expected_device_name, isolation)
                if run.get("reported_devices") and "signature" in run:
                    label = run["reported_devices"][0]
                    if device_labels.setdefault(settings["backend"], label) != label:
                        run["error"] = f"GPU Device changed during the benchmark: {label}"
                signature = run.get("signature")
                if repeat == 1 and name == "baseline":
                    group["reference"] = signature if "error" not in run else None
                baseline = group["reference"]
                run["valid"] = "error" not in run and baseline is not None and signature == baseline
                if not run["valid"]:
                    run.setdefault("error", "PDBQT or stdout pose table differs from the same-card baseline")
                if a100_signatures:
                    run["a100_comparison"] = compare_signatures(signature, a100_signatures[case_id])
                group["runs"][name].append(run)
                print(f"{repeat}/{args.repeats} {case_id} {name} {run['wall_s']:.3f}s {run.get('error', 'completed')}",
                      flush=True)

    summarize(report, versions)
    result = output / "result.json"
    result.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"{result} (valid={report['valid']})", flush=True)
    return 0 if report["valid"] else 1


def main():
    parser = argparse.ArgumentParser(description="Compare two Vina-GPU builds on one Moore Threads GPU (Linux, Python 3.10+)")
    parser.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    parser.add_argument("--baseline", type=Path, required=True, help="Local baseline executable")
    parser.add_argument("--candidate", type=Path, required=True, help="Local candidate executable")
    parser.add_argument("--backend", choices=("opencl", "musa"), required=True, help="Baseline backend")
    parser.add_argument("--candidate-backend", choices=("opencl", "musa"), help="Defaults to --backend")
    parser.add_argument("--baseline-kernels", type=Path, help="OpenCL kernel binaries (defaults to baseline executable directory)")
    parser.add_argument("--candidate-kernels", type=Path, help="OpenCL kernel binaries (defaults to candidate executable directory)")
    parser.add_argument("--gpu", required=True, help="Verified single GPU identity, preferably model, PCI address and driver")
    parser.add_argument("--device-probe", help="Optional command, without shell expansion, that prints one line per visible GPU")
    parser.add_argument("--expected-device-name", help="Expected GPU Device label in the stdout of both versions")
    parser.add_argument("--cpu-affinity", help="Comma-separated Linux CPU IDs (defaults to inherited affinity)")
    parser.add_argument("--cache-env", default=[], action="append", metavar="NAME",
                        help="Set a runtime-specific cache-path environment variable independently for each version")
    parser.add_argument("--isolate", action="store_true", help="Run each case in a network-isolated bwrap sandbox")
    parser.add_argument("--bwrap", help="Absolute bwrap executable path (default /usr/bin/bwrap)")
    parser.add_argument("--device-node", action="append", default=[], help="GPU character device to expose with --isolate")
    parser.add_argument("--library-dir", action="append", default=[], help="Runtime lib directory to expose with --isolate")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--reference", type=Path, help="A100 result JSON for cross-hardware output comparison")
    parser.add_argument("--output", type=Path, required=True, help="New results directory")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if not args.gpu.strip():
        parser.error("--gpu must identify the verified single GPU")
    if args.isolate and not args.device_node:
        parser.error("--isolate requires an explicit --device-node")
    if not args.isolate and (args.bwrap or args.device_node or args.library_dir):
        parser.error("--bwrap, --device-node and --library-dir require --isolate")
    for name in args.cache_env:
        if not re.fullmatch(r"[A-Z][A-Z0-9_]*", name) or name in {
            "HOME", "LC_ALL", "TMPDIR", "XDG_CACHE_HOME", *DEVICE_ENV, *STREAM_ENV
        }:
            parser.error(f"Invalid cache path variable: {name}")
    try:
        with tempfile.TemporaryDirectory(dir="/dev/shm", prefix="vina-astra-") as tmp:
            return benchmark(args, Path(tmp))
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"Benchmark setup failed: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
