#!/usr/bin/env python3
"""
EdgeFEM Simulation: test_wg_sim
Type: waveguide

Run with: python run_simulation.py
"""

import sys
from pathlib import Path
import numpy as np

# Add EdgeFEM to path
edgefem_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(edgefem_root / "build" / "python"))

import pyedgefem as em


def main():
    # Configuration
    geo_file = "waveguide.geo"
    msh_file = geo_file.replace(".geo", ".msh")

    # Mesh the geometry (if not already done)
    if not Path(msh_file).exists():
        import subprocess
        print(f"Meshing {geo_file}...")
        result = subprocess.run(["gmsh", geo_file, "-3", "-o", msh_file],
                                capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Meshing failed: {result.stderr}")
            return

    # Load mesh
    print(f"Loading mesh from {msh_file}...")
    mesh = em.load_gmsh(msh_file)
    print(f"  Nodes: {len(mesh.nodes)}")
    print(f"  Tets: {len(mesh.tets)}")
    print(f"  Tris: {len(mesh.tris)}")

    # Validate mesh
    quality = em.validate_mesh(mesh)
    if quality.is_valid():
        print("  Mesh quality: OK")
    else:
        print(f"  Warning: {quality.num_degenerate_tets} degenerate, "
              f"{quality.num_inverted_tets} inverted tets")

    # Set up boundary conditions
    # Adjust PEC tag based on your geometry
    bc = em.build_edge_pec(mesh, 1)
    print(f"  Applied PEC boundary")

    # Frequency setup
    freq = 10e9  # Modify as needed
    omega = 2 * np.pi * freq

    params = em.MaxwellParams()
    params.omega = omega
    params.eps_r = 1.0 + 0j
    params.mu_r = 1.0 + 0j

    # Solve
    print(f"\nSolving at f = {freq/1e9:.2f} GHz...")

    # For simulations without ports, just assemble and solve
    asm = em.assemble_maxwell(mesh, params, bc, [], -1)

    opts = em.SolveOptions()
    opts.tolerance = 1e-10
    opts.max_iterations = 10000

    result = em.solve_linear(asm.A, asm.b, opts)

    print(f"  Converged: {result.converged}")
    print(f"  Iterations: {result.iters}")
    print(f"  Residual: {result.residual:.2e}")

    # TODO: Add your post-processing here
    # - S-parameter extraction
    # - Far-field patterns
    # - Field visualization

    print("\nSimulation complete!")


if __name__ == "__main__":
    main()
