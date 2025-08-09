#!/usr/bin/env bash
set -euo pipefail
gmsh -3 rect_waveguide.geo -format msh2 -o rect_waveguide.msh
