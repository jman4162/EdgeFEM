#!/usr/bin/env bash
set -euo pipefail
build_dir=${1:-build}
cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" -j
"$build_dir/src/vectorem_scalar_demo" examples/cube_cavity.msh 1
