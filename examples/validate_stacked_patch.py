#!/usr/bin/env python3
"""
Stacked Patch Antenna Validation Script

Validates the end-to-end simulation pipeline: mesh generation (OCC kernel),
lumped port excitation, Maxwell assembly, solve, and S-parameter extraction.

Simulates a single-patch cavity-fed antenna at 2.45 GHz on FR-4 substrate
and checks for physically reasonable results.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "python"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from edgefem.designs.stacked_patch import StackedPatchDesign

C0 = 299792458.0
MU0 = 4e-7 * np.pi


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


def run_validation():
    """Run single-patch validation."""
    print("=" * 70)
    print("EdgeFEM Stacked Patch Antenna Validation")
    print("=" * 70)

    W = 38e-3
    L = 29e-3
    h = 1.6e-3
    eps_r = 4.4
    tan_d = 0.02
    cavity_size = 45e-3
    cavity_depth = 5e-3
    probe_radius = 5e-3  # 5mm for mesh resolution at density=10
    probe_y_offset = -L / 6
    Z0 = 50.0

    print(f"\nDesign Parameters:")
    print(f"  Patch: {W*1e3:.1f} x {L*1e3:.1f} mm")
    print(f"  Substrate: {h*1e3:.1f} mm FR-4 (eps_r={eps_r}, tan_d={tan_d})")
    print(f"  Cavity: {cavity_size*1e3:.1f} mm, {cavity_depth*1e3:.1f} mm deep")
    print(f"  Probe: radius={probe_radius*1e3:.1f} mm, offset=(0, {probe_y_offset*1e3:.2f}) mm")

    f_r_ana, eps_eff, L_eff = analytical_resonant_freq(L, W, h, eps_r)
    print(f"\nAnalytical: f_r = {f_r_ana/1e9:.4f} GHz, eps_eff = {eps_eff:.3f}")

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
    print("\nGenerating mesh (density=10)...")
    try:
        design.generate_mesh(density=10, design_freq=2.45e9)
    except RuntimeError as e:
        print(f"ERROR: Mesh generation failed: {e}")
        print("Ensure gmsh Python module is installed: pip install gmsh")
        return False

    print(f"  Mesh: {design.mesh.num_nodes()} nodes, "
          f"{design.mesh.num_tets()} tets, {design.mesh.num_edges()} edges")

    # Single frequency simulation
    results = {}

    print("\nSimulating at design frequency...")
    try:
        S11, _ = design.simulate(2.45e9)
        s11_mag = abs(S11)
        s11_dB = 20 * np.log10(max(s11_mag, 1e-30))
        Z_in = Z0 * (1 + S11) / (1 - S11) if abs(1 - S11) > 1e-30 else complex(1e6, 0)

        print(f"  |S11| = {s11_mag:.4f} ({s11_dB:.1f} dB)")
        print(f"  Z_in = {Z_in.real:.1f} + j{Z_in.imag:.1f} ohms")

        results["solver_converged"] = True
        results["passivity"] = s11_mag <= 1.01
        results["s11_below_unity"] = s11_mag < 1.01
    except Exception as e:
        print(f"  ERROR: {e}")
        results["solver_converged"] = False
        return False

    # Frequency sweep (5 points for speed)
    print("\nFrequency sweep: 2.0 - 3.0 GHz (5 points)...")
    try:
        freqs, Z_in_arr, S11_arr = design.frequency_sweep(
            2.0e9, 3.0e9, 5, verbose=True
        )
        results["sweep_completed"] = True

        print(f"\n{'Freq (GHz)':>12} {'|S11| (dB)':>12} {'R_in':>12} {'X_in':>12}")
        print("-" * 52)
        for i in range(len(freqs)):
            s_dB = 20 * np.log10(max(abs(S11_arr[i]), 1e-30))
            print(f"  {freqs[i]/1e9:10.4f} {s_dB:12.2f} "
                  f"{Z_in_arr[i].real:12.1f} {Z_in_arr[i].imag:12.1f}")

        # Passivity check
        max_s11 = np.max(np.abs(S11_arr))
        results["all_passive"] = max_s11 <= 1.01

        # Check for S11 dip (resonance) — at least one point below -1 dB
        min_s11_dB = 20 * np.log10(max(np.min(np.abs(S11_arr)), 1e-30))
        results["has_resonance"] = min_s11_dB < -1.0
        print(f"\n  Min |S11| = {min_s11_dB:.1f} dB")

        # Check impedance at best match point
        best_idx = np.argmin(np.abs(S11_arr))
        R_in = Z_in_arr[best_idx].real
        X_in = Z_in_arr[best_idx].imag
        results["reasonable_impedance"] = (10 < R_in < 300 and abs(X_in) < 200)
        print(f"  Best match at {freqs[best_idx]/1e9:.3f} GHz: "
              f"Z_in = {R_in:.1f} + j{X_in:.1f} ohms")
    except Exception as e:
        print(f"  ERROR: {e}")
        results["sweep_completed"] = False

    # Summary
    print("\n" + "=" * 70)
    print("Validation Summary")
    print("=" * 70)

    checks = [
        ("Mesh generation", True),
        ("Solver convergence", results.get("solver_converged", False)),
        ("Passivity (|S11| <= 1)", results.get("passivity", False)),
        ("S11 < 1 (port coupling)", results.get("s11_below_unity", False)),
        ("Frequency sweep", results.get("sweep_completed", False)),
        ("Resonance dip (< -1 dB)", results.get("has_resonance", False)),
        ("Reasonable impedance", results.get("reasonable_impedance", False)),
    ]

    all_pass = True
    for name, passed in checks:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  [{status}] {name}")

    if all_pass:
        print("\nAll validation checks passed.")
    else:
        print("\nSome checks failed.")

    return all_pass


if __name__ == "__main__":
    success = run_validation()
    sys.exit(0 if success else 1)
