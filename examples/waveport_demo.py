#!/usr/bin/env python3
"""
WavePort Demo - End-to-end example of wave port S-parameter extraction.

This example demonstrates the complete workflow for computing S-parameters
of a rectangular waveguide using modal wave ports.

Requirements:
    - VectorEM built with VECTOREM_PYTHON=ON
    - Mesh file: examples/rect_waveguide.msh
"""

import numpy as np

try:
    import pyvectorem as em
except ImportError:
    print("Error: pyvectorem not found. Build VectorEM with -DVECTOREM_PYTHON=ON")
    raise


def main():
    """Run complete wave port S-parameter extraction workflow."""

    # --------------------------------------------------------------------------
    # Step 1: Load the mesh
    # --------------------------------------------------------------------------
    mesh_path = "examples/rect_waveguide.msh"
    print(f"Loading mesh from {mesh_path}...")
    mesh = em.load_gmsh(mesh_path)
    print(f"  Loaded {len(mesh.nodes)} nodes, {len(mesh.tets)} tets")

    # --------------------------------------------------------------------------
    # Step 2: Define boundary conditions
    # --------------------------------------------------------------------------
    # Physical tag 1 = PEC walls (perfect electric conductor)
    pec_tag = 1
    bc = em.build_edge_pec(mesh, pec_tag)
    print(f"  Applied PEC boundary on physical tag {pec_tag}")

    # --------------------------------------------------------------------------
    # Step 3: Set up wave ports
    # --------------------------------------------------------------------------
    # For a WR-90 waveguide: a=22.86mm, b=10.16mm
    # TE10 mode cutoff: fc = c/(2a) = 6.56 GHz

    # Port surface tags (from Gmsh geometry)
    port1_tag = 2  # Input port
    port2_tag = 3  # Output port

    # Operating frequency (must be above cutoff for TE10)
    freq = 10e9  # 10 GHz
    omega = 2 * np.pi * freq

    # Solve for TE10 mode on each port
    print(f"\nSolving TE10 mode at f = {freq/1e9:.2f} GHz...")

    # Extract port surface meshes
    port1_surface = em.extract_surface_mesh(mesh, port1_tag)
    port2_surface = em.extract_surface_mesh(mesh, port2_tag)

    # Create TE10 mode for WR-90 waveguide
    # Modal properties are computed analytically for rectangular waveguides
    eps_r = 1.0 + 0j
    mu_r = 1.0 + 0j

    # Solve 2D eigenmode on port cross-section
    # (In practice, solve_port_eigens extracts modes from the surface mesh)
    num_modes = 1
    modes1 = em.solve_port_eigens(port1_surface.mesh, num_modes, omega, eps_r, mu_r,
                                   em.ModePolarization.TE)
    modes2 = em.solve_port_eigens(port2_surface.mesh, num_modes, omega, eps_r, mu_r,
                                   em.ModePolarization.TE)

    if len(modes1) == 0 or len(modes2) == 0:
        print("Error: No modes found. Check mesh and frequency.")
        return

    mode1 = modes1[0]
    mode2 = modes2[0]
    print(f"  Port 1: fc = {mode1.fc/1e9:.3f} GHz, Z0 = {mode1.Z0.real:.2f} ohms")
    print(f"  Port 2: fc = {mode2.fc/1e9:.3f} GHz, Z0 = {mode2.Z0.real:.2f} ohms")

    # Build wave port structures (project modal fields onto 3D edges)
    wave_port1 = em.build_wave_port(mesh, port1_surface, mode1)
    wave_port2 = em.build_wave_port(mesh, port2_surface, mode2)
    ports = [wave_port1, wave_port2]

    print(f"  Port 1: {len(wave_port1.edges)} edges")
    print(f"  Port 2: {len(wave_port2.edges)} edges")

    # --------------------------------------------------------------------------
    # Step 4: Set up Maxwell parameters
    # --------------------------------------------------------------------------
    params = em.MaxwellParams()
    params.omega = omega
    params.eps_r = eps_r
    params.mu_r = mu_r

    # --------------------------------------------------------------------------
    # Step 5: Calculate S-parameters
    # --------------------------------------------------------------------------
    print("\nCalculating S-parameters...")
    S = em.calculate_sparams(mesh, params, bc, ports)

    print(f"\nS-parameter matrix ({S.shape[0]}x{S.shape[1]}):")
    print(f"  S11 = {S[0,0]:.4f} ({20*np.log10(abs(S[0,0])):.2f} dB)")
    print(f"  S21 = {S[1,0]:.4f} ({20*np.log10(abs(S[1,0])):.2f} dB)")
    print(f"  S12 = {S[0,1]:.4f} ({20*np.log10(abs(S[0,1])):.2f} dB)")
    print(f"  S22 = {S[1,1]:.4f} ({20*np.log10(abs(S[1,1])):.2f} dB)")

    # Check passivity
    power_in = np.abs(S[0,0])**2 + np.abs(S[1,0])**2
    print(f"\n  Passivity check: |S11|^2 + |S21|^2 = {power_in:.6f}")
    if power_in <= 1.0 + 1e-6:
        print("  PASS: Network is passive")
    else:
        print("  WARN: Network appears non-passive (numerical error?)")

    # --------------------------------------------------------------------------
    # Step 6: Frequency sweep (optional)
    # --------------------------------------------------------------------------
    print("\nRunning frequency sweep...")
    freqs = np.linspace(8e9, 12e9, 11)
    s11_data = []
    s21_data = []

    for f in freqs:
        params.omega = 2 * np.pi * f

        # Re-solve modes at each frequency for accurate dispersion
        modes1 = em.solve_port_eigens(port1_surface.mesh, 1, params.omega,
                                       eps_r, mu_r, em.ModePolarization.TE)
        modes2 = em.solve_port_eigens(port2_surface.mesh, 1, params.omega,
                                       eps_r, mu_r, em.ModePolarization.TE)

        if len(modes1) > 0 and len(modes2) > 0:
            wp1 = em.build_wave_port(mesh, port1_surface, modes1[0])
            wp2 = em.build_wave_port(mesh, port2_surface, modes2[0])
            S = em.calculate_sparams(mesh, params, bc, [wp1, wp2])
            s11_data.append((f, S[0,0]))
            s21_data.append((f, S[1,0]))
            print(f"  f = {f/1e9:.2f} GHz: S11 = {20*np.log10(abs(S[0,0])):.2f} dB, "
                  f"S21 = {20*np.log10(abs(S[1,0])):.2f} dB")

    # --------------------------------------------------------------------------
    # Step 7: Export results
    # --------------------------------------------------------------------------
    # Write to Touchstone format
    output_path = "examples/waveport_demo_sparams.s2p"
    sparams_list = []
    for i, f in enumerate(freqs):
        if i < len(s11_data):
            s2 = em.SParams2()
            s2.s11 = s11_data[i][1]
            s2.s21 = s21_data[i][1]
            s2.s12 = s21_data[i][1]  # Reciprocal
            s2.s22 = s11_data[i][1]  # Symmetric
            sparams_list.append(s2)

    if sparams_list:
        em.write_touchstone(output_path, list(freqs[:len(sparams_list)]), sparams_list)
        print(f"\nS-parameters written to {output_path}")

    print("\nDemo complete.")


if __name__ == "__main__":
    main()
