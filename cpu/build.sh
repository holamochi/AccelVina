#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/rebuilds"
dst=$(mktemp -d "$here/rebuilds/build-XXXXXXXX")
mkdir "$dst/tmp"
export TMPDIR="$dst/tmp"
ln -s "$here/src" "$dst/src"
make -C "$dst" -f "$here/Makefile" -j "${JOBS:-4}" BUILD=build-02 PROFILE_DIR="$here/profiles" all 2>&1 | tee "$dst/build.log"
if rg -Fq 'profile count data file not found' "$dst/build.log"; then
    echo 'GCC did not read all measured profile files' >&2
    exit 1
fi
sha256sum "$dst/build-02/vina"
