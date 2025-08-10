#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

gmsh -3 cube_cavity.geo -format msh2 -o cube_cavity.msh
gmsh -3 mixed_cavity.geo -format msh2 -o mixed_cavity.msh

