# Examples

This directory contains simple meshes and scripts for smoke testing.

## cube_cavity

Generate the mesh (requires [Gmsh](https://gmsh.info/)):

```bash
./make_mesh.sh
```

This creates `cube_cavity.msh` with exterior triangle faces tagged as physical `1`.

## rect_waveguide

Simulate a straight WR-90 rectangular waveguide and sweep S-parameters.

1. Generate the mesh (requires [Gmsh](https://gmsh.info/)):

   ```bash
   ./make_waveguide_mesh.sh
   ```

   View the mesh with `gmsh rect_waveguide.msh` if desired.

2. Build and run the demo sweep:

   ```bash
   ./run_waveguide_demo.sh
   ```

   This writes `examples/waveguide_sparams.csv`.
   
3. Plot, export Touchstone, and compare against the analytic solution:

   ```bash
   python examples/plot_waveguide_sparams.py
   ```

   The script generates a plot image (`examples/waveguide_sparams.png`,
   not tracked in git) and writes a Touchstone file
   `examples/waveguide_sparams.s2p`.

Sample CSV and Touchstone outputs are included in this repository for
reference; regenerate them with the commands above to verify the
results locally.

## patch_antenna_pattern

`patch_antenna_pattern.csv` shows a sample 2D far-field pattern exported
from a patch antenna simulation. Visualize it with the helper script:

```bash
python ../tools/plot_pattern.py patch_antenna_pattern.csv
```

This produces a polar plot of the gain pattern in the $xz$-plane.

## patch_antenna_design

`patch_antenna_design.py` performs a simple analytic design of a
2.45 GHz inset‑fed rectangular microstrip patch on FR‑4. It sweeps the
inset depth and patch length, exporting S‑parameter data and the final
geometry dimensions.

Run:

```bash
python patch_antenna_design.py
```

The script writes `patch_sparams.csv`, a Touchstone file
`patch_2p45_FR4.s1p`, a `patch_antenna.geo` model, and a
`patch_design.json` summary in this directory.

To generate a coarse tetrahedral mesh and run the nascent full-wave solver,
use:

```bash
PYTHONPATH=../build/python ./run_patch_fullwave.py
```

This writes `patch_antenna.msh` and a placeholder residual sweep
`patch_fullwave_sparams.csv`.
