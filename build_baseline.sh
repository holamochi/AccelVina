#!/usr/bin/env bash
set -euo pipefail
ulimit -c 0

if [[ $# != 2 ]]; then
    echo "Usage: bash build_baseline.sh NEW_BUILD_DIR DEPENDENCY_ARCHIVE_DIR" >&2
    exit 2
fi

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build=$(realpath -m -- "$1")
archives=$(realpath -- "$2")
source_dir="$root/upstream/vina-gpu-2.1"
example_dir="$root/benchmarks/examples/gpu-2bm2"
commit=180272b8a5265d6ed9664178345933cebe2cd349
(cd -- "$source_dir" && sha256sum --quiet -c SHA256SUMS)

(cd -- "$archives" && sha256sum -c <<'HASHES'
e5eeb98b824cb3a680e09f1e2c4d9274493361da37fd561bf51071ca917253fa  libboost-1.84.0-h6fcfa73_0.conda
f6a6eb40a33b32eaab5d9ab36567b126f54c0b2112fe53131cf7b55afff2d0d9  libboost-headers-1.84.0-ha770c72_0.conda
159f2a550592bae49859fee83d372acd152328fdf95c0dcd8b9409f8fad5db93  OpenCL-Headers-v2024.10.24.tar.gz
HASHES
)

mkdir -- "$build"
mkdir -- "$build/src" "$build/deps" "$build/runtime" "$build/runtime/lib" "$build/runtime/icd"
cp -a -- "$source_dir"/{Makefile,lib,main,OpenCL} "$build/src/"
python - "$archives" "$build/deps" <<'PY'
import sys
from pathlib import Path

from conda_package_handling.api import extract

archives, deps = map(Path, sys.argv[1:])
extract(str(archives / "libboost-1.84.0-h6fcfa73_0.conda"), dest_dir=str(deps / "boost-libs"))
extract(str(archives / "libboost-headers-1.84.0-ha770c72_0.conda"), dest_dir=str(deps / "boost-headers"))
PY
tar -xzf "$archives/OpenCL-Headers-v2024.10.24.tar.gz" -C "$build/deps"

runtime="$build/runtime"
for library in program_options system filesystem thread atomic chrono; do
    cp -L -- "$build/deps/boost-libs/lib/libboost_${library}.so.1.84.0" "$runtime/lib/"
done
cp -L -- /usr/lib/x86_64-linux-gnu/libOpenCL.so.1 "$runtime/lib/"
printf '%s\n' libnvidia-opencl.so.1 > "$runtime/icd/nvidia.icd"

# The container also has CUDA forward-compat libraries from a newer driver.
export LD_LIBRARY_PATH="$runtime/lib:/usr/lib/x86_64-linux-gnu"
export OCL_ICD_VENDORS="$runtime/icd"
export CUDA_CACHE_DISABLE=1
cd -- "$build/src"
make_args=(
    "WORK_DIR=$PWD"
    "BOOST_INC_PATH=-I$build/deps/boost-headers/include"
    "OPENCL_INC_PATH=-I$build/deps/OpenCL-Headers-2024.10.24"
    "LIB_PATH=-L$runtime/lib"
    'LIB1=-l:libboost_program_options.so.1.84.0 -l:libboost_system.so.1.84.0 -l:libboost_filesystem.so.1.84.0 -l:libboost_thread.so.1.84.0 -l:libOpenCL.so.1'
    'SRC=./lib/*.cpp ./OpenCL/src/wrapcl.cpp'
    'OPENCL_VERSION=-DOPENCL_3_0'
    'GPU_PLATFORM=-DNVIDIA_PLATFORM'
    'DOCKING_BOX_SIZE=-DSMALL_BOX'
)
make source "${make_args[@]}" 2>&1 | tee "$build/build-kernel.log"
./AutoDock-Vina-GPU-2-1 \
    --receptor "$example_dir/2bm2_protein.pdbqt" \
    --ligand "$example_dir/2bm2_ligand.pdbqt" \
    --center_x 40.415 --center_y 110.986 --center_z 82.673 \
    --size_x 30 --size_y 30 --size_z 30 --thread 1000 --seed 20261012 \
    --out "$build/kernel-smoke.pdbqt" --opencl_binary_path . \
    > "$build/kernel-smoke.log" 2>&1
test -s "$build/kernel-smoke.pdbqt"
test -s Kernel1_Opt.bin && test -s Kernel2_Opt.bin
cp -- Kernel1_Opt.bin Kernel2_Opt.bin "$runtime/"

curl -fsSL --retry 2 \
    "https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU-2.1/$commit/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1" \
    -o "$runtime/AutoDock-Vina-GPU-2-1"
(cd -- "$runtime" && sha256sum -c <<'HASHES'
9d27a4e88f3dd597297168d74943e7f8a151d7986ffb118f6b30b9bd349826b1  AutoDock-Vina-GPU-2-1
HASHES
)
chmod +x "$runtime/AutoDock-Vina-GPU-2-1"
{
    printf 'upstream_commit=%s\n' "$commit"
    printf 'executable=upstream tracked Linux binary\n'
    readelf -p .comment "$runtime/AutoDock-Vina-GPU-2-1"
    printf 'kernel generation host compiler\n'
    gcc --version | head -n 1
    nvidia-smi --query-gpu=uuid,name,driver_version --format=csv,noheader
} > "$build/build-environment.txt"
(cd -- "$runtime" && sha256sum AutoDock-Vina-GPU-2-1 Kernel{1,2}_Opt.bin lib/*.so* icd/nvidia.icd > SHA256SUMS)
