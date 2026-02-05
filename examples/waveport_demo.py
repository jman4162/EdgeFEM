#!/usr/bin/env python3
"""
WavePort Demo - S-parameter extraction workflow for rectangular waveguide.

This example demonstrates the VectorEM workflow for computing S-parameters:
1. Load mesh and apply boundary conditions
2. Extract port surfaces and compute modal fields
3. Build wave port structures
4. Calculate S-parameters (analytical or full-wave FEM)
5. Export to Touchstone format

Requirements:
    - VectorEM built with VECTOREM_PYTHON=ON
    - Mesh file: examples/rect_waveguide.msh

Generate mesh with:
    gmsh examples/rect_waveguide.geo -3 -format msh2 -o examples/rect_waveguide.msh
"""

import numpy as np
from pathlib import Path

try:
    import pyvectorem as em
except ImportError:
    print("Error: pyvectorem not found. Build VectorEM with -DVECTOREM_PYTHON=ON")
    raise


def main():
    """Run wave port S-parameter extraction workflow."""

    script_dir = Path(__file__).parent

    # ==========================================================================
    # STEP 1: Load Mesh
    # ==========================================================================
    mesh_path = script_dir / "rect_waveguide.msh"
    print(f"Loading mesh from {mesh_path}...")
    mesh = em.load_gmsh(str(mesh_path))
    print(f"  Nodes: {len(mesh.nodes)}, Tets: {len(mesh.tets)}")

    # Validate mesh quality
    quality = em.validate_mesh(mesh)
    if quality.is_valid():
        print("  Mesh quality: OK")
    else:
        print(f"  Warning: {quality.num_degenerate_tets} degenerate, "
              f"{quality.num_inverted_tets} inverted tets")

    # ==========================================================================
    # STEP 2: Apply Boundary Conditions
    # ==========================================================================
    # Physical groups from Gmsh:
    #   Tag 1 = PEC walls
    #   Tag 2 = Port1 (z=0)
    #   Tag 3 = Port2 (z=L)
    pec_tag = 1
    bc = em.build_edge_pec(mesh, pec_tag)
    print(f"  Applied PEC boundary on tag {pec_tag}")

    # ==========================================================================
    # STEP 3: Extract Port Surfaces
    # ==========================================================================
    port1_surface = em.extract_surface_mesh(mesh, 2)  # Input port
    port2_surface = em.extract_surface_mesh(mesh, 3)  # Output port
    print(f"  Port surfaces: {len(port1_surface.mesh.tris)} triangles each")

    # ==========================================================================
    # STEP 4: Define Waveguide Parameters
    # ==========================================================================
    # WR-90 rectangular waveguide
    a = 0.02286  # width = 22.86 mm
    b = 0.01016  # height = 10.16 mm
    L = 0.05     # length = 50 mm
    port = em.RectWaveguidePort(a, b)

    # TE10 cutoff frequency: fc = c/(2a) ≈ 6.56 GHz
    c0 = 299792458
    fc_te10 = c0 / (2 * a)
    print(f"\nWR-90 Waveguide: a = {a*1000:.2f} mm, b = {b*1000:.2f} mm")
    print(f"  TE10 cutoff: fc = {fc_te10/1e9:.3f} GHz")

    # ==========================================================================
    # STEP 5: Compute Analytical S-Parameters
    # ==========================================================================
    # For a straight, lossless, matched waveguide section:
    # - S11 = S22 = 0 (no reflection)
    # - S21 = S12 = exp(-j*beta*L) (transmission with phase delay)

    print("\n" + "="*60)
    print("ANALYTICAL S-PARAMETERS (Straight Matched Waveguide)")
    print("="*60)

    freqs = np.linspace(7e9, 12e9, 11)
    S_matrices = []

    for f in freqs:
        # Get analytical mode properties
        mode = em.solve_te10_mode(port, f)

        # Get analytical S-parameters
        s = em.straight_waveguide_sparams(port, L, f)

        # Build full S-matrix
        S = np.array([[s.s11, s.s12], [s.s21, s.s22]], dtype=complex)
        S_matrices.append(S)

        if f < fc_te10:
            status = "EVANESCENT"
        else:
            status = f"beta = {mode.beta.real:6.2f} rad/m"

        s21_db = 20*np.log10(max(abs(s.s21), 1e-15))
        print(f"  f = {f/1e9:5.2f} GHz: |S21| = {abs(s.s21):.4f} ({s21_db:6.2f} dB)  {status}")

    # ==========================================================================
    # STEP 6: Export to Touchstone Format
    # ==========================================================================
    output_path = script_dir / "waveport_analytical.s2p"
    em.write_touchstone_nport(str(output_path), list(freqs), S_matrices)
    print(f"\nWritten: {output_path}")

    # ==========================================================================
    # STEP 7: Demonstrate Coupling Matrix Conversion
    # ==========================================================================
    print("\n" + "="*60)
    print("COUPLING MATRIX (S → Z → Y)")
    print("="*60)

    # At 10 GHz (center of passband)
    f_center = 10e9
    mode = em.solve_te10_mode(port, f_center)
    s = em.straight_waveguide_sparams(port, L, f_center)
    S = np.array([[s.s11, s.s12], [s.s21, s.s22]], dtype=complex)

    # Convert to coupling matrix (Z0 = modal impedance)
    coupling = em.compute_coupling(S, f_center, mode.Z0.real)

    print(f"\nAt f = {f_center/1e9:.1f} GHz (Z0 = {mode.Z0.real:.1f} ohms):")
    print(f"\n  S-matrix:")
    print(f"    S11 = {S[0,0]: .4f}    S12 = {S[0,1]: .4f}")
    print(f"    S21 = {S[1,0]: .4f}    S22 = {S[1,1]: .4f}")

    print(f"\n  Z-matrix (ohms):")
    print(f"    Z11 = {coupling.Z[0,0]: .1f}    Z12 = {coupling.Z[0,1]: .1f}")
    print(f"    Z21 = {coupling.Z[1,0]: .1f}    Z22 = {coupling.Z[1,1]: .1f}")

    print(f"\n  Y-matrix (mS):")
    print(f"    Y11 = {coupling.Y[0,0]*1000: .3f}    Y12 = {coupling.Y[0,1]*1000: .3f}")
    print(f"    Y21 = {coupling.Y[1,0]*1000: .3f}    Y22 = {coupling.Y[1,1]*1000: .3f}")

    # Passivity check
    is_pass = em.is_passive(S)
    print(f"\n  Passivity check: {'PASS' if is_pass else 'FAIL'}")

    # ==========================================================================
    # STEP 8: Show Full-Wave FEM Workflow (Reference)
    # ==========================================================================
    print("\n" + "="*60)
    print("FULL-WAVE FEM WORKFLOW (Reference)")
    print("="*60)
    print("""
The full-wave FEM workflow follows these steps:

1. Compute analytical mode:
   mode = em.solve_te10_mode(port, freq)

2. Populate mode field on mesh:
   em.populate_te10_field(surface, port, mode)

3. Build wave port structure:
   wave_port = em.build_wave_port(mesh, surface, mode)

4. Set Maxwell parameters:
   params = em.MaxwellParams()
   params.omega = 2 * pi * freq
   params.eps_r = 1.0 + 0j
   params.mu_r = 1.0 + 0j

5. Calculate S-parameters:
   S = em.calculate_sparams(mesh, params, bc, [port1, port2])

Note: Full-wave FEM results should match analytical for simple
waveguide geometries. Discrepancies indicate mesh or solver
issues that need investigation.
""")

    print("Demo complete.")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
