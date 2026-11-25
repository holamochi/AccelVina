from __future__ import annotations

import os
from pathlib import Path


KERNELS = ("Kernel1_Opt.bin", "Kernel2_Opt.bin")
STAGED_LIGANDS_ROOT = Path("/dev/shm")
RUNTIME_ENV = (
    "LC_ALL", "LD_LIBRARY_PATH", "LD_PRELOAD", "MUSA_VISIBLE_DEVICES", "CUDA_VISIBLE_DEVICES",
    "OCL_ICD_VENDORS", "OCL_ICD_FILENAMES", "OPENCL_VENDOR_PATH", "GPU_DEVICE_ORDINAL",
    "VINA_EXPECT_PCI", "CUDA_CACHE_DISABLE", "CUDA_CACHE_MAXSIZE", "VINA_MUSA_PROFILE",
)


def _absolute(path):
    path = Path(path)
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Expected an absolute path without parent components: {path}")
    return path


def _musa_lib_mount():
    sdk = Path("/usr/local/musa")
    lib = sdk / "lib"
    if not lib.is_dir():
        return None

    # Preserve the link in the read-only /usr mount. Only the lib subtree of an
    # SDK installed outside /usr is needed for the executable's SDK RPATH.
    for path, suffix in ((sdk, "lib"), (sdk.resolve(strict=True) / "lib", "")):
        while path.is_symlink():
            path = Path(os.path.normpath(path.parent / os.readlink(path)))
            if not path.is_relative_to("/usr"):
                return lib.resolve(strict=True), path / suffix
    return None


def isolated_command(command, work, env, *, device_nodes, library_dirs=(), cache_env=(),
                     bwrap="/usr/bin/bwrap"):
    """Build the argv for one frozen case. Run it with a clean outer Popen env."""
    if len(command) < 3 or len(command) % 2 != 1:
        raise ValueError("Expected the executable followed by benchmark option/value pairs")
    options = dict(zip(command[1::2], command[2::2], strict=True))
    if len(options) != (len(command) - 1) // 2 or not all(key.startswith("--") for key in options):
        raise ValueError("Invalid or duplicate benchmark option")
    if "--receptor" not in options or ("--ligand" in options) == ("--ligand_directory" in options):
        raise ValueError("Expected a receptor and exactly one ligand input")

    work = _absolute(work)
    argv = [str(_absolute(bwrap)), "--unshare-all", "--die-with-parent", "--new-session",
            "--cap-drop", "ALL", "--clearenv", "--proc", "/proc", "--dev", "/dev",
            "--tmpfs", "/dev/shm", "--tmpfs", "/tmp", "--ro-bind", "/usr", "/usr",
            "--ro-bind", "/sys", "/sys"]

    def bind(path, *, directory=False, writable=False, target=None):
        path = _absolute(path)
        source = path.resolve(strict=True)
        if not (source.is_dir() if directory else source.is_file()):
            raise ValueError(f"Unexpected bind source type: {source}")
        argv.extend(("--bind" if writable else "--ro-bind", str(source), str(target or path)))

    for path in map(Path, ("/lib", "/lib64", "/bin")):
        if path.is_symlink():
            argv.extend(("--symlink", os.readlink(path), str(path)))
        elif path.is_dir():
            bind(path, directory=True)
    for path in map(Path, ("/etc/ld.so.cache", "/etc/nsswitch.conf", "/etc/passwd", "/etc/group")):
        if path.is_file():
            bind(path)

    if sdk_mount := _musa_lib_mount():
        source, target = sdk_mount
        bind(source, directory=True, target=target)
    for path in map(_absolute, library_dirs):
        if path.name not in {"lib", "lib64"} or path.resolve(strict=True).name not in {"lib", "lib64"}:
            raise ValueError(f"Bind the specific runtime library directory: {path}")
        bind(path, directory=True)
    for path in map(_absolute, device_nodes):
        if not path.is_relative_to("/dev") or not path.is_char_device():
            raise ValueError(f"Expected an explicit GPU character device under /dev: {path}")
        argv.extend(("--dev-bind", str(path.resolve(strict=True)), str(path)))

    bind(command[0])
    for option in ("--receptor", "--flex", "--ligand"):
        if option in options:
            bind(options[option])
    if "--ligand_directory" in options:
        ordered = _absolute(options["--ligand_directory"])
        staging_root = STAGED_LIGANDS_ROOT.resolve(strict=True)
        resolved = ordered.resolve(strict=True)
        if (ordered.is_symlink() or not ordered.is_dir() or resolved == staging_root
                or not resolved.is_relative_to(staging_root)):
            raise ValueError(f"Expected the runner's staged ligand directory: {ordered}")
        entries = list(ordered.iterdir())
        if not entries or any(p.is_symlink() or not p.is_file() or p.suffix != ".pdbqt" for p in entries):
            raise ValueError(f"Expected only staged PDBQT files in {ordered}")
        bind(ordered, directory=True)
    if "--opencl_binary_path" in options:
        for name in KERNELS:
            bind(Path(options["--opencl_binary_path"]) / name)
        icd_dir = Path("/etc/OpenCL/vendors")
        if icd_dir.is_dir():
            bind(icd_dir, directory=True)

    writable = {work / "out"}
    cache_root = work.parents[2] / "driver-cache" / work.name
    child_env = {name: env[name] for name in RUNTIME_ENV if name in env}
    for name in ("HOME", "XDG_CACHE_HOME", *cache_env):
        value = env.get(name)
        if value is None:
            raise ValueError(f"Missing runner cache directory for {name}")
        path = _absolute(value)
        if path == cache_root or not path.is_relative_to(cache_root):
            raise ValueError(f"Cache directory outside this version's scope: {path}")
        writable.add(path)
        child_env[name] = value

    for path in sorted(writable):
        if path.resolve(strict=True) != path:
            raise ValueError(f"Writable mount cannot be a symlink: {path}")
        bind(path, directory=True, writable=True)

    if profile := child_env.get("VINA_MUSA_PROFILE"):
        if not any(_absolute(profile).is_relative_to(path) for path in writable):
            raise ValueError("VINA_MUSA_PROFILE must be within an isolated writable directory")
    child_env.update(PATH="/usr/bin:/bin", TMPDIR="/tmp")
    for name, value in sorted(child_env.items()):
        argv.extend(("--setenv", name, value))
    return [*argv, "--chdir", str(work), "--", *command]
