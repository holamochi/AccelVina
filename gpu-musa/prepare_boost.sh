#!/usr/bin/env bash
# Extract the two upstream-pinned conda-forge artifacts into a private prefix.
set -euo pipefail
[[ $# == 1 || $# == 2 ]] || { echo "usage: prepare_boost.sh NEW_PREFIX [ARCHIVE_DIR]" >&2; exit 2; }
prefix=$(realpath -m "$1")
mkdir "$prefix"
archives=${2:-"$prefix/archives"}
mkdir -p "$archives"
archives=$(realpath "$archives")
cd "$archives"
verify_extract() {
    local name=$1 digest=$2
    if [[ ! -f $name ]]; then
        curl -fsSL --retry 2 "https://conda.anaconda.org/conda-forge/linux-64/$name" -o "$name"
    fi
    echo "$digest  $name" | sha256sum -c -
    unzip -p "$name" 'pkg-*.tar.zst' | tar --zstd -xf - -C "$prefix"
}
verify_extract libboost-1.84.0-h6fcfa73_0.conda e5eeb98b824cb3a680e09f1e2c4d9274493361da37fd561bf51071ca917253fa
verify_extract libboost-headers-1.84.0-ha770c72_0.conda f6a6eb40a33b32eaab5d9ab36567b126f54c0b2112fe53131cf7b55afff2d0d9
