#!/usr/bin/env python3
"""
Example 02: Frequency Sweep Analysis

This example demonstrates broadband S-parameter analysis across the WR-90
waveguide passband (7-12 GHz), computing S-parameters at multiple frequencies
and plotting the results.

Steps:
1. Load mesh and set up ports
2. Sweep frequency from 7 to 12 GHz
3. Compute S-parameters at each frequency point
4. Plot |S11| and |S21| vs frequency
5. Export to Touchstone format

Requirements:
    - VectorEM built with VECTOREM_PYTHON=ON
    - matplotlib for plotting (optional)
    - Mesh file: examples/rect_waveguide.msh
"""

import numpy as np
from pathlib import Path
import sys

try:
    import pyvectorem as em
except ImportError:
    print("Error: pyvectorem not found. Build VectorEM with -DVECTOREM_PYTHON=ON")
    raise


def main():
    """Run frequency sweep analysis."""

    script_dir = Path(__file__).parent

    print("=" * 60)
    print("VectorEM Example 02: Frequency Sweep Analysis")
    print("=" * 60)

    # ========================================================================
    # STEP 1: Load Mesh and Set Up Ports
    # ========================================================================
    print("\n1. Setting up simulation...")

    mesh_path = script_dir / "rect_waveguide.msh"
    mesh = em.load_gmsh(str(mesh_path))
    bc = em.build_edge_pec(mesh, 1)

    # WR-90 waveguide
    a = 0.02286  # width
    b = 0.01016  # height
    L = 0.05     # length
    port = em.RectWaveguidePort(a, b)

    # Extract port surfaces
    port1_surface = em.extract_surface_mesh(mesh, 2)
    port2_surface = em.extract_surface_mesh(mesh, 3)

    # Compute eigenvector (shared across frequencies for same mode)
    kc_te10 = np.pi / a
    kc_sq = kc_te10 ** 2
    pec_edges_set = set(bc.dirichlet_edges)
    eigenvector = em.compute_te_eigenvector(mesh, pec_edges_set, kc_sq)

    print(f"   Mesh: {len(mesh.tets)} tets, {len(mesh.edges)} edges")
    print(f"   Eigenvector computed")

    # ========================================================================
    # STEP 2: Define Frequency Sweep
    # ========================================================================
    print("\n2. Defining frequency sweep...")

    c0 = 299792458
    fc_te10 = c0 / (2 * a)

    # Sweep from 7 to 12 GHz (11 points)
    freq_start = 7e9
    freq_stop = 12e9
    n_points = 11
    frequencies = np.linspace(freq_start, freq_stop, n_points)

    print(f"   TE10 cutoff: {fc_te10/1e9:.2f} GHz")
    print(f"   Frequency range: {freq_start/1e9:.1f} - {freq_stop/1e9:.1f} GHz")
    print(f"   Number of points: {n_points}")

    # ========================================================================
    # STEP 3: Compute S-Parameters at Each Frequency
    # ========================================================================
    print("\n3. Computing S-parameters...")
    print(f"   {'Freq(GHz)':<12} {'|S11|':<10} {'|S21|':<10} {'Phase(S21)':<12} Status")
    print("   " + "-" * 52)

    S_matrices = []
    s11_data = []
    s21_data = []
    s21_phase_data = []

    for freq in frequencies:
        # Build ports for this frequency
        mode1 = em.solve_te10_mode(port, freq)
        mode2 = em.solve_te10_mode(port, freq)

        wp1 = em.build_wave_port_from_eigenvector(mesh, port1_surface, eigenvector, mode1, pec_edges_set)
        wp2 = em.build_wave_port_from_eigenvector(mesh, port2_surface, eigenvector, mode2, pec_edges_set)

        # Set up parameters
        params = em.MaxwellParams()
        params.omega = 2 * np.pi * freq
        params.port_abc_scale = 0.5

        # Calculate S-parameters
        S = em.calculate_sparams_eigenmode(mesh, params, bc, [wp1, wp2])
        S_matrices.append(S)

        # Extract values
        s11 = abs(S[0, 0])
        s21 = abs(S[1, 0])
        phase_s21 = np.angle(S[1, 0]) * 180 / np.pi

        s11_data.append(s11)
        s21_data.append(s21)
        s21_phase_data.append(phase_s21)

        # Status
        status = "OK" if s21 > 0.9 and s11 < 0.15 else "CHECK"
        print(f"   {freq/1e9:<12.2f} {s11:<10.4f} {s21:<10.4f} {phase_s21:<12.1f} {status}")

    # ========================================================================
    # STEP 4: Plot Results (if matplotlib available)
    # ========================================================================
    print("\n4. Generating plots...")

    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        freq_ghz = frequencies / 1e9

        # |S11| in dB
        ax1 = axes[0, 0]
        s11_db = 20 * np.log10(np.maximum(s11_data, 1e-6))
        ax1.plot(freq_ghz, s11_db, 'b-o', linewidth=2, markersize=6)
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('|S11| (dB)')
        ax1.set_title('Return Loss')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([-40, 0])

        # |S21| in dB
        ax2 = axes[0, 1]
        s21_db = 20 * np.log10(np.maximum(s21_data, 1e-6))
        ax2.plot(freq_ghz, s21_db, 'r-o', linewidth=2, markersize=6)
        ax2.set_xlabel('Frequency (GHz)')
        ax2.set_ylabel('|S21| (dB)')
        ax2.set_title('Insertion Loss')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim([-1, 0.1])

        # |S21| linear
        ax3 = axes[1, 0]
        ax3.plot(freq_ghz, s21_data, 'g-o', linewidth=2, markersize=6)
        ax3.axhline(y=0.95, color='k', linestyle='--', alpha=0.5, label='Target (0.95)')
        ax3.set_xlabel('Frequency (GHz)')
        ax3.set_ylabel('|S21|')
        ax3.set_title('Transmission Magnitude')
        ax3.grid(True, alpha=0.3)
        ax3.set_ylim([0.9, 1.01])
        ax3.legend()

        # Phase
        ax4 = axes[1, 1]
        ax4.plot(freq_ghz, s21_phase_data, 'm-o', linewidth=2, markersize=6)
        ax4.set_xlabel('Frequency (GHz)')
        ax4.set_ylabel('Phase(S21) (degrees)')
        ax4.set_title('Transmission Phase')
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_path = script_dir / "02_frequency_sweep.png"
        plt.savefig(plot_path, dpi=150)
        print(f"   Plot saved: {plot_path}")
        plt.close()

    except ImportError:
        print("   matplotlib not available, skipping plots")

    # ========================================================================
    # STEP 5: Export to Touchstone Format
    # ========================================================================
    print("\n5. Exporting to Touchstone format...")
    output_path = script_dir / "02_frequency_sweep.s2p"
    em.write_touchstone_nport(str(output_path), list(frequencies), S_matrices)
    print(f"   Written: {output_path}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Frequency range: {freq_start/1e9:.1f} - {freq_stop/1e9:.1f} GHz")
    print(f"Number of points: {n_points}")
    print(f"Average |S21|: {np.mean(s21_data):.4f}")
    print(f"Min |S21|: {np.min(s21_data):.4f} at {frequencies[np.argmin(s21_data)]/1e9:.1f} GHz")
    print(f"Max |S11|: {np.max(s11_data):.4f} at {frequencies[np.argmax(s11_data)]/1e9:.1f} GHz")

    all_pass = all(s21 > 0.9 and s11 < 0.15 for s11, s21 in zip(s11_data, s21_data))
    print(f"\nAll points meet criteria (|S21| > 0.9, |S11| < 0.15): {'YES' if all_pass else 'NO'}")

    print("\nDemo complete.")


if __name__ == "__main__":
    main()
