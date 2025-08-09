#!/usr/bin/env bash
set -euo pipefail
build_dir=${1:-build}

if command -v clang-format >/dev/null; then
  clang-format --dry-run --Werror $(git ls-files '*.cpp' '*.hpp')
else
  echo "clang-format not found; skipping format check" >&2
fi

if command -v clang-tidy >/dev/null; then
  clang-tidy -p "$build_dir" -warnings-as-errors='*' $(git ls-files 'src/*.cpp')
else
  echo "clang-tidy not found; skipping analysis" >&2
fi

