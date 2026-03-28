#!/usr/bin/env python3
"""
Stacked Patch Antenna Case Study

Generates validation plots and a summary report for the cavity-fed
stacked patch antenna simulation using EdgeFEM lumped port excitation.

Outputs:
  examples/stacked_patch_s11.png      - |S11| vs frequency
  examples/stacked_patch_impedance.png - R_in and X_in vs frequency
  examples/stacked_patch_smith.png    - Smith chart
  examples/stacked_patch_report.txt   - Text summary

Usage:
  PYTHONPATH=build/python:python python3 examples/stacked_patch_case_study.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "python"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from edgefem.designs.stacked_patch import StackedPatchDesign

C0 = 299792458.0
MU0 = 4e-7 * np.pi
OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def analytical_resonant_freq(L, W, h, eps_r):
    """Transmission-line model resonant frequency (Balanis 14.2)."""
    eps_eff = (eps_r + 1) / 2 + (eps_r - 1) / 2 * (1 + 12 * h / W) ** (-0.5)
    delta_L = 0.412 * h * (
        (eps_eff + 0.3) * (W / h + 0.264)
        / ((eps_eff - 0.258) * (W / h + 0.8))
    )
    L_eff = L + 2 * delta_L
    f_r = C0 / (2 * L_eff * np.sqrt(eps_eff))
    return f_r, eps_eff, L_eff


def plot_s11(freqs_ghz, S11_arr, f_r_ana_ghz, save_path):
    """Plot |S11| in dB vs frequency."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s11_dB = 20 * np.log10(np.maximum(np.abs(S11_arr), 1e-30))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(freqs_ghz, s11_dB, "b-o", linewidth=2, markersize=4, label="EdgeFEM")
    ax.axvline(f_r_ana_ghz, color="r", linestyle="--", alpha=0.7,
               label=f"Analytical f_r = {f_r_ana_ghz:.3f} GHz")
    ax.axhline(-10, color="gray", linestyle=":", alpha=0.5, label="-10 dB")

    # Mark minimum
    idx_min = np.argmin(s11_dB)
    ax.plot(freqs_ghz[idx_min], s11_dB[idx_min], "rv", markersize=10,
            label=f"Min: {s11_dB[idx_min]:.1f} dB @ {freqs_ghz[idx_min]:.3f} GHz")

    ax.set_xlabel("Frequency (GHz)", fontsize=12)
    ax.set_ylabel("|S11| (dB)", fontsize=12)
    ax.set_title("Stacked Patch Antenna — Return Loss", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(min(s11_dB) - 3, 3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {save_path}")


def plot_impedance(freqs_ghz, Z_in_arr, save_path):
    """Plot R_in and X_in vs frequency."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(freqs_ghz, Z_in_arr.real, "b-o", linewidth=2, markersize=4, label="R_in")
    ax.plot(freqs_ghz, Z_in_arr.imag, "r-s", linewidth=2, markersize=4, label="X_in")
    ax.axhline(50, color="gray", linestyle=":", alpha=0.5, label="50 ohm")
    ax.axhline(0, color="black", linestyle="-", alpha=0.2)

    ax.set_xlabel("Frequency (GHz)", fontsize=12)
    ax.set_ylabel("Impedance (ohm)", fontsize=12)
    ax.set_title("Stacked Patch Antenna — Input Impedance", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {save_path}")


def plot_smith(S11_arr, save_path):
    """Plot S11 on a Smith chart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 7))

    # Draw Smith chart circles
    theta = np.linspace(0, 2 * np.pi, 200)

    # Unit circle (|Gamma|=1)
    ax.plot(np.cos(theta), np.sin(theta), "k-", linewidth=1.5)

    # Constant resistance circles: r = 0, 0.5, 1, 2
    for r in [0, 0.5, 1, 2, 5]:
        center_x = r / (1 + r)
        radius = 1 / (1 + r)
        ax.plot(center_x + radius * np.cos(theta),
                radius * np.sin(theta), "gray", linewidth=0.5, alpha=0.5)

    # Constant reactance arcs: x = +/-0.5, +/-1, +/-2
    for x in [0.5, 1, 2]:
        center_y = 1 / x
        radius = 1 / x
        arc = center_y + radius * np.exp(1j * theta)
        mask = np.abs(arc) <= 1.01
        real_pts = arc.real
        imag_pts = arc.imag
        for sign in [1, -1]:
            cx, cy = 1.0, sign * center_y
            pts_x = cx + radius * np.cos(theta)
            pts_y = cy + radius * np.sin(theta)
            inside = pts_x**2 + pts_y**2 <= 1.02
            pts_x[~inside] = np.nan
            pts_y[~inside] = np.nan
            ax.plot(pts_x, pts_y, "gray", linewidth=0.5, alpha=0.5)

    # Horizontal axis
    ax.plot([-1, 1], [0, 0], "gray", linewidth=0.5)

    # Plot S11 data
    ax.plot(S11_arr.real, S11_arr.imag, "b-o", linewidth=2, markersize=5)
    ax.plot(S11_arr.real[0], S11_arr.imag[0], "go", markersize=8, label="Start")
    ax.plot(S11_arr.real[-1], S11_arr.imag[-1], "rs", markersize=8, label="End")

    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect("equal")
    ax.set_title("Stacked Patch Antenna — Smith Chart", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {save_path}")


def generate_report(params, freqs, Z_in_arr, S11_arr, f_r_ana, mesh_info, save_path):
    """Generate text summary report."""
    s11_dB = 20 * np.log10(np.maximum(np.abs(S11_arr), 1e-30))
    idx_min = np.argmin(s11_dB)
    f_res_fem = freqs[idx_min]

    lines = []
    lines.append("=" * 70)
    lines.append("EdgeFEM Stacked Patch Antenna Case Study")
    lines.append("=" * 70)
    lines.append("")
    lines.append("1. DESIGN PARAMETERS")
    lines.append("-" * 40)
    lines.append(f"  Patch width:       {params['W']*1e3:.1f} mm")
    lines.append(f"  Patch length:      {params['L']*1e3:.1f} mm")
    lines.append(f"  Substrate height:  {params['h']*1e3:.1f} mm")
    lines.append(f"  Substrate eps_r:   {params['eps_r']}")
    lines.append(f"  Loss tangent:      {params['tan_d']}")
    lines.append(f"  Cavity size:       {params['cavity_size']*1e3:.1f} mm")
    lines.append(f"  Cavity depth:      {params['cavity_depth']*1e3:.1f} mm")
    lines.append(f"  Probe radius:      {params['probe_radius']*1e3:.1f} mm")
    lines.append(f"  Probe Y offset:    {params['probe_y_offset']*1e3:.2f} mm")
    lines.append(f"  Reference Z0:      {params['Z0']} ohm")
    lines.append("")
    lines.append("2. MESH")
    lines.append("-" * 40)
    lines.append(f"  Nodes:             {mesh_info['nodes']}")
    lines.append(f"  Tetrahedra:        {mesh_info['tets']}")
    lines.append(f"  Edges (DOFs):      {mesh_info['edges']}")
    lines.append(f"  Mesh density:      {params['density']} elements/wavelength")
    lines.append("")
    lines.append("3. ANALYTICAL REFERENCE")
    lines.append("-" * 40)
    lines.append(f"  f_r (TL model):    {f_r_ana/1e9:.4f} GHz")
    lines.append("")
    lines.append("4. FEM RESULTS")
    lines.append("-" * 40)
    lines.append(f"  Frequency range:   {freqs[0]/1e9:.2f} - {freqs[-1]/1e9:.2f} GHz")
    lines.append(f"  Number of points:  {len(freqs)}")
    lines.append(f"  Best match freq:   {f_res_fem/1e9:.4f} GHz")
    lines.append(f"  Min |S11|:         {s11_dB[idx_min]:.1f} dB")
    lines.append(f"  Z_in at resonance: {Z_in_arr[idx_min].real:.1f} + "
                 f"j{Z_in_arr[idx_min].imag:.1f} ohm")
    lines.append(f"  Max |S11|:         {np.max(np.abs(S11_arr)):.4f} (passivity)")
    lines.append("")
    lines.append("5. FREQUENCY SWEEP DATA")
    lines.append("-" * 40)
    lines.append(f"  {'Freq (GHz)':>12} {'|S11| (dB)':>12} {'R_in':>10} {'X_in':>10}")
    for i in range(len(freqs)):
        lines.append(f"  {freqs[i]/1e9:12.4f} {s11_dB[i]:12.2f} "
                     f"{Z_in_arr[i].real:10.1f} {Z_in_arr[i].imag:10.1f}")
    lines.append("")
    lines.append("6. COMPARISON: ANALYTICAL vs FEM")
    lines.append("-" * 40)
    lines.append(f"  Analytical f_r:    {f_r_ana/1e9:.4f} GHz")
    lines.append(f"  FEM f_r:           {f_res_fem/1e9:.4f} GHz")
    shift_pct = (f_res_fem - f_r_ana) / f_r_ana * 100
    lines.append(f"  Shift:             {shift_pct:+.1f}%")
    lines.append("")
    lines.append("  Note: The FEM resonance may shift from the analytical TL model")
    lines.append("  prediction due to: (1) cavity-backing lowers the effective")
    lines.append("  resonant frequency, (2) coarse mesh (~10 elem/lambda) affects")
    lines.append("  accuracy, and (3) frequency-independent normalization is an")
    lines.append("  approximation. The cavity effect typically dominates, pulling")
    lines.append("  the resonance lower than the free-space TL model predicts.")
    lines.append("")
    lines.append("7. KNOWN LIMITATIONS")
    lines.append("-" * 40)
    lines.append("  - Coarse mesh (10 elem/lambda) shifts resonance upward")
    lines.append("  - Impedance normalization uses frequency-independent factor")
    lines.append("  - Radiation patterns not available (Huygens surface not meshed)")
    lines.append("  - Single probe feed; dual-polarization not modeled")
    lines.append("")
    lines.append("8. VALIDATION CHECKS")
    lines.append("-" * 40)
    all_passive = np.all(np.abs(S11_arr) <= 1.01)
    has_resonance = np.min(s11_dB) < -1.0
    reasonable_Z = (10 < Z_in_arr[idx_min].real < 300 and
                    abs(Z_in_arr[idx_min].imag) < 200)
    checks = [
        ("Passivity (|S11| <= 1)", all_passive),
        ("Resonance dip (< -1 dB)", has_resonance),
        ("Reasonable impedance", reasonable_Z),
    ]
    for name, passed in checks:
        lines.append(f"  [{'PASS' if passed else 'FAIL'}] {name}")
    lines.append("")
    lines.append("=" * 70)

    text = "\n".join(lines)
    with open(save_path, "w") as f:
        f.write(text)
    print(f"  Saved: {save_path}")
    print()
    print(text)


def main():
    print("EdgeFEM Stacked Patch Antenna Case Study")
    print("=" * 50)

    # Design parameters
    W = 38e-3
    L = 29e-3
    h = 1.6e-3
    eps_r = 4.4
    tan_d = 0.02
    cavity_size = 45e-3
    cavity_depth = 5e-3
    probe_radius = 5e-3
    probe_y_offset = -L / 6
    Z0 = 50.0
    density = 10

    params = dict(W=W, L=L, h=h, eps_r=eps_r, tan_d=tan_d,
                  cavity_size=cavity_size, cavity_depth=cavity_depth,
                  probe_radius=probe_radius, probe_y_offset=probe_y_offset,
                  Z0=Z0, density=density)

    # Analytical reference
    f_r_ana, eps_eff, L_eff = analytical_resonant_freq(L, W, h, eps_r)
    print(f"\nAnalytical f_r = {f_r_ana/1e9:.4f} GHz (eps_eff={eps_eff:.3f})")

    # Create design
    design = StackedPatchDesign(
        patches=[{"width": W, "length": L, "height": h}],
        substrates=[{"thickness": h, "eps_r": eps_r, "tan_d": tan_d}],
        cavity_size=cavity_size,
        cavity_depth=cavity_depth,
        probe_radius=probe_radius,
        probe_offset=(0.0, probe_y_offset),
        z0=Z0,
    )

    # Generate mesh
    print(f"\nGenerating mesh (density={density})...")
    design.generate_mesh(density=density, design_freq=2.45e9)
    mesh_info = {
        "nodes": design.mesh.num_nodes(),
        "tets": design.mesh.num_tets(),
        "edges": design.mesh.num_edges(),
    }
    print(f"  {mesh_info['nodes']} nodes, {mesh_info['tets']} tets, "
          f"{mesh_info['edges']} edges")

    # Frequency sweep
    n_points = 11
    f_start, f_stop = 2.0e9, 3.5e9
    print(f"\nFrequency sweep: {f_start/1e9:.1f} - {f_stop/1e9:.1f} GHz "
          f"({n_points} points)...")
    freqs, Z_in_arr, S11_arr = design.frequency_sweep(
        f_start, f_stop, n_points, verbose=True
    )

    # Generate plots
    print("\nGenerating plots...")
    plot_s11(
        freqs / 1e9, S11_arr, f_r_ana / 1e9,
        os.path.join(OUT_DIR, "stacked_patch_s11.png"),
    )
    plot_impedance(
        freqs / 1e9, Z_in_arr,
        os.path.join(OUT_DIR, "stacked_patch_impedance.png"),
    )
    plot_smith(
        S11_arr,
        os.path.join(OUT_DIR, "stacked_patch_smith.png"),
    )

    # Export Touchstone
    s1p_path = os.path.join(OUT_DIR, "stacked_patch.s1p")
    design.export_touchstone(s1p_path, freqs, S11_arr)
    print(f"  Saved: {s1p_path}")

    # Generate report
    print("\nGenerating report...")
    generate_report(
        params, freqs, Z_in_arr, S11_arr, f_r_ana, mesh_info,
        os.path.join(OUT_DIR, "stacked_patch_report.txt"),
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
