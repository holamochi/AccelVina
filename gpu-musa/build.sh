#!/usr/bin/env bash
# build.sh baseline|candidate NEW_DIR BOOST_PREFIX [UPSTREAM_SOURCE]
set -euo pipefail
[[ $# -ge 3 && $# -le 4 && ($1 == baseline || $1 == candidate) ]] || {
    echo "usage: build.sh baseline|candidate NEW_DIR BOOST_PREFIX [UPSTREAM_SOURCE]" >&2; exit 2; }
project=$(cd "$(dirname "$0")" && pwd)
kind=$1
out=$(realpath -m "$2")
boost=$(realpath "$3")
src=$(realpath "${4:-$project/../Vina-GPU-2.1/AutoDock-Vina-GPU-2.1}")
sdk=$(realpath "${MUSA_HOME:-/usr/local/musa}")
mkdir "$out"
mkdir "$out/src" "$out/lib"
for name in program_options filesystem thread atomic; do
    cp -L "$boost/lib/libboost_${name}.so.1.84.0" "$out/lib/"
done
cp -a "$src/lib" "$src/main" "$src/OpenCL" "$out/src/"
cp -a "$project/port" "$project/runtime" "$project/build.mk" "$out/"
if [[ $kind == candidate ]]; then
    mkdir "$out/reference"
    cp "$src/OpenCL/src/kernels/kernel1.cl" "$out/reference/"
    cp "$project/port/kernel1.mu" "$out/reference/grid_fallback.mu"
    patch --batch --fuzz=0 -d "$out/src" -p1 < "$project/opt/0001-host-data-preparation.patch"
    for delta in "$project"/opt/*-musa-*.patch; do
        patch --batch --fuzz=0 -d "$out" -p1 < "$delta"
    done
fi
cd "$out"
make -f build.mk -j"${JOBS:-4}" BOOST="$boost" SDK="$sdk" ARCH="${MUSA_ARCH:-mp_31}"
sha256sum AutoDock-Vina-GPU-2-1 obj/kernel*.o lib/*.so.* > SHA256SUMS
