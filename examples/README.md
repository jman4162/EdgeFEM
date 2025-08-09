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

3. Plot and compare against the analytic solution:

   ```bash
   python examples/plot_waveguide_sparams.py
   ```

   The plot is saved as `examples/waveguide_sparams.png`.

