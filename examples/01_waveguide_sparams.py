#!/usr/bin/env python3
"""
Example 01: Waveguide S-Parameter Extraction

This example demonstrates the basic EdgeFEM workflow for computing S-parameters
of a WR-90 rectangular waveguide using the eigenmode-based method.

Steps:
1. Load mesh and apply boundary conditions
2. Extract port surfaces
3. Compute 3D FEM eigenvector for TE10 mode
4. Build wave ports with eigenvector-based weights
5. Calculate S-parameters using eigenmode method
6. Export to Touchstone format

Requirements:
    - EdgeFEM built with EDGEFEM_PYTHON=ON
    - Mesh file: examples/rect_waveguide.msh

Generate mesh with:
    gmsh examples/rect_waveguide.geo -3 -format msh2 -o examples/rect_waveguide.msh
"""

import numpy as np
from pathlib import Path

try:
    import pyedgefem as em
except ImportError:
    print("Error: pyedgefem not found. Build EdgeFEM with -DEDGEFEM_PYTHON=ON")
    raise


def main():
    """Run waveguide S-parameter extraction workflow."""

    script_dir = Path(__file__).parent

    print("=" * 60)
    print("EdgeFEM Example 01: Waveguide S-Parameter Extraction")
    print("=" * 60)

    # ========================================================================
    # STEP 1: Load Mesh
    # ========================================================================
    mesh_path = script_dir / "rect_waveguide.msh"
    print(f"\n1. Loading mesh from {mesh_path}...")
    mesh = em.load_gmsh(str(mesh_path))
    print(f"   Nodes: {len(mesh.nodes)}")
    print(f"   Tetrahedra: {len(mesh.tets)}")
    print(f"   Edges: {len(mesh.edges)}")

    # ========================================================================
    # STEP 2: Apply Boundary Conditions
    # ========================================================================
    print("\n2. Applying PEC boundary condition on waveguide walls...")
    # Physical groups from Gmsh:
    #   Tag 1 = PEC walls
    #   Tag 2 = Port1 (input)
    #   Tag 3 = Port2 (output)
    pec_tag = 1
    bc = em.build_edge_pec(mesh, pec_tag)
    print(f"   PEC edges: {len(bc.dirichlet_edges)}")

    # ========================================================================
    # STEP 3: Define Waveguide Parameters
    # ========================================================================
    print("\n3. Defining WR-90 waveguide parameters...")
    # WR-90 rectangular waveguide dimensions
    a = 0.02286  # width = 22.86 mm
    b = 0.01016  # height = 10.16 mm
    L = 0.05     # length = 50 mm
    port = em.RectWaveguidePort(a, b)

    # Physical constants
    c0 = 299792458  # speed of light (m/s)
    fc_te10 = c0 / (2 * a)  # TE10 cutoff frequency

    print(f"   Width: {a * 1000:.2f} mm")
    print(f"   Height: {b * 1000:.2f} mm")
    print(f"   Length: {L * 1000:.2f} mm")
    print(f"   TE10 cutoff: {fc_te10 / 1e9:.3f} GHz")

    # ========================================================================
    # STEP 4: Extract Port Surfaces
    # ========================================================================
    print("\n4. Extracting port surfaces...")
    port1_surface = em.extract_surface_mesh(mesh, 2)  # Input port
    port2_surface = em.extract_surface_mesh(mesh, 3)  # Output port
    print(f"   Port 1: {len(port1_surface.mesh.tris)} triangles")
    print(f"   Port 2: {len(port2_surface.mesh.tris)} triangles")

    # ========================================================================
    # STEP 5: Compute 3D FEM Eigenvector for TE10 Mode
    # ========================================================================
    print("\n5. Computing 3D FEM eigenvector for TE10 mode...")
    kc_te10 = np.pi / a  # Cutoff wavenumber
    kc_sq = kc_te10 ** 2

    # Get PEC edges as a set for eigenvector computation
    pec_edges_set = set(bc.dirichlet_edges)
    eigenvector = em.compute_te_eigenvector(mesh, pec_edges_set, kc_sq)
    print(f"   Eigenvector computed (norm = {np.linalg.norm(eigenvector):.4f})")

    # ========================================================================
    # STEP 6: Build Wave Ports with Eigenvector-Based Weights
    # ========================================================================
    print("\n6. Building wave ports with eigenvector-based weights...")
    freq = 10e9  # Operating frequency: 10 GHz

    # Compute analytical mode parameters
    mode1 = em.solve_te10_mode(port, freq)
    mode2 = em.solve_te10_mode(port, freq)

    print(f"   Frequency: {freq / 1e9:.1f} GHz")
    print(f"   TE10 impedance: {mode1.Z0.real:.1f} ohms")
    print(f"   Propagation constant: {mode1.beta.real:.2f} rad/m")

    # Build ports using eigenvector
    wp1 = em.build_wave_port_from_eigenvector(mesh, port1_surface, eigenvector, mode1, pec_edges_set)
    wp2 = em.build_wave_port_from_eigenvector(mesh, port2_surface, eigenvector, mode2, pec_edges_set)
    print(f"   Port 1: {len(wp1.edges)} edges")
    print(f"   Port 2: {len(wp2.edges)} edges")

    # ========================================================================
    # STEP 7: Calculate S-Parameters
    # ========================================================================
    print("\n7. Calculating S-parameters...")

    # Set up Maxwell parameters with optimal ABC scaling
    params = em.MaxwellParams()
    params.omega = 2 * np.pi * freq
    params.port_abc_scale = 0.5  # Optimal value for accurate S-parameters

    # Calculate using eigenmode method
    S = em.calculate_sparams_eigenmode(mesh, params, bc, [wp1, wp2])

    print("\n   S-Parameter Matrix:")
    print(f"   S11 = {S[0,0]:.4f}  (|S11| = {abs(S[0,0]):.4f})")
    print(f"   S21 = {S[1,0]:.4f}  (|S21| = {abs(S[1,0]):.4f})")
    print(f"   S12 = {S[0,1]:.4f}  (|S12| = {abs(S[0,1]):.4f})")
    print(f"   S22 = {S[1,1]:.4f}  (|S22| = {abs(S[1,1]):.4f})")

    # Expected values
    beta = mode1.beta.real
    expected_phase = -beta * L * 180 / np.pi
    while expected_phase < -180:
        expected_phase += 360
    while expected_phase > 180:
        expected_phase -= 360

    actual_phase = np.angle(S[1, 0]) * 180 / np.pi

    print(f"\n   Expected: |S11| = 0, |S21| = 1, phase(S21) = {expected_phase:.1f}deg")
    print(f"   Actual:   |S11| = {abs(S[0,0]):.3f}, |S21| = {abs(S[1,0]):.3f}, phase(S21) = {actual_phase:.1f}deg")

    # Passivity check
    passivity = abs(S[0, 0])**2 + abs(S[1, 0])**2
    print(f"\n   Passivity check: |S11|^2 + |S21|^2 = {passivity:.4f} (should be <= 1)")

    # ========================================================================
    # STEP 8: Export to Touchstone Format
    # ========================================================================
    print("\n8. Exporting to Touchstone format...")
    output_path = script_dir / "01_waveguide_sparams.s2p"

    # For single frequency, create a simple list
    freqs = [freq]
    S_matrices = [S]

    em.write_touchstone_nport(str(output_path), freqs, S_matrices)
    print(f"   Written: {output_path}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Waveguide: WR-90, {L*1000:.0f}mm length")
    print(f"Frequency: {freq/1e9:.1f} GHz")
    print(f"Results:")
    print(f"  - Transmission: |S21| = {abs(S[1,0])*100:.1f}%")
    print(f"  - Reflection: |S11| = {abs(S[0,0])*100:.1f}%")
    print(f"  - Phase error: {abs(actual_phase - expected_phase):.1f} degrees")
    print(f"  - Passivity: {'PASS' if passivity <= 1.01 else 'FAIL'}")

    print("\nDemo complete.")


if __name__ == "__main__":
    main()
