#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
build_dir=${1:-build}
cmake -S "$script_dir/.." -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" -j
"$build_dir/src/vectorem_waveguide_demo" > "$script_dir/waveguide_sparams.csv"
