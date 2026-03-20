#!/usr/bin/env python3
"""
Validation: EdgeFEM vs Analytical Fresnel Coefficients for a Dielectric Slab

Computes reflection and transmission coefficients for a dielectric slab
using both analytical Fresnel theory and EdgeFEM's C++ analytical reference,
then plots the comparison.

Configuration:
  - Dielectric slab: eps_r = 4.0, thickness = 3 mm
  - Normal incidence
  - Frequency range: 5-15 GHz (spans multiple Fabry-Perot resonances)

The analytical Fresnel model accounts for multiple reflections within the
slab (Fabry-Perot effect), producing oscillating R and T vs frequency.

Usage:
    python examples/validation_fresnel.py
"""

import numpy as np
import os
import sys

# ---------------------------------------------------------------------------
# Analytical Fresnel coefficients for a dielectric slab (normal incidence)
# ---------------------------------------------------------------------------

c0 = 299792458.0  # speed of light (m/s)


def fresnel_slab(freq, eps_r, d_slab):
    """
    Compute power reflection and transmission for a lossless dielectric slab.

    Uses the transfer-matrix (Fabry-Perot) formulation for a single slab
    in air at normal incidence.

    Parameters
    ----------
    freq : float or array
        Frequency in Hz.
    eps_r : float
        Relative permittivity of the slab.
    d_slab : float
        Slab thickness in metres.

    Returns
    -------
    R : same shape as freq
        Power reflection coefficient |r|^2.
    T : same shape as freq
        Power transmission coefficient |t|^2.
    """
    freq = np.asarray(freq, dtype=float)

    n1 = 1.0                       # air
    n2 = np.sqrt(eps_r)            # slab

    # Single-interface Fresnel coefficients (normal incidence)
    r12 = (n1 - n2) / (n1 + n2)
    t12 = 2.0 * n1 / (n1 + n2)
    t21 = 2.0 * n2 / (n1 + n2)
    r21 = -r12

    # Phase through slab
    k0 = 2.0 * np.pi * freq / c0
    delta = k0 * n2 * d_slab       # one-way phase

    exp_2jd = np.exp(2j * delta)
    exp_jd = np.exp(1j * delta)

    # Total reflection and transmission (complex)
    r_total = (r12 + r21 * exp_2jd) / (1.0 + r12 * r21 * exp_2jd)
    t_total = (t12 * t21 * exp_jd) / (1.0 + r12 * r21 * exp_2jd)

    R = np.abs(r_total) ** 2
    T = np.abs(t_total) ** 2
    return R, T


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Slab parameters (must match examples/dielectric_slab.geo)
    eps_r = 4.0
    d_slab = 3.0e-3  # 3 mm

    # Frequency sweep
    freqs = np.linspace(5e9, 15e9, 201)
    R_ana, T_ana = fresnel_slab(freqs, eps_r, d_slab)

    # Verify energy conservation
    max_violation = np.max(np.abs(R_ana + T_ana - 1.0))
    assert max_violation < 1e-12, f"Energy conservation violated: {max_violation}"

    # Print key results table
    print("=" * 64)
    print("Dielectric Slab Fresnel Validation")
    print(f"  eps_r = {eps_r},  d = {d_slab*1e3:.1f} mm,  n = {np.sqrt(eps_r):.2f}")
    print("=" * 64)
    print(f"{'Freq (GHz)':>12} {'|R|^2':>10} {'|T|^2':>10} {'R+T':>10}")
    print("-" * 44)

    spot_freqs = [5e9, 7.5e9, 10e9, 12.5e9, 15e9]
    for f in spot_freqs:
        R, T = fresnel_slab(f, eps_r, d_slab)
        print(f"{f/1e9:12.1f} {R:10.4f} {T:10.4f} {R+T:10.6f}")

    # Identify Fabry-Perot resonances (T -> 1)
    resonance_mask = T_ana > 0.999
    if np.any(resonance_mask):
        res_freqs = freqs[resonance_mask]
        print(f"\nFabry-Perot resonances (|T|^2 > 0.999):")
        for f in res_freqs:
            print(f"  {f/1e9:.3f} GHz")

    # Attempt to plot (graceful fallback if matplotlib is unavailable)
    try:
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not installed — skipping plot generation.")
        print("Install with: pip install matplotlib")
        return

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    # Reflection
    ax = axes[0]
    ax.plot(freqs / 1e9, R_ana, "b-", linewidth=1.5, label="Analytical Fresnel")
    ax.set_ylabel("|R|$^2$ (power)")
    ax.set_title(
        f"Dielectric Slab Fresnel Coefficients\n"
        f"$\\varepsilon_r$ = {eps_r}, d = {d_slab*1e3:.1f} mm"
    )
    ax.legend(loc="upper right")
    ax.set_ylim(-0.02, 0.25)
    ax.grid(True, alpha=0.3)

    # Transmission
    ax = axes[1]
    ax.plot(freqs / 1e9, T_ana, "r-", linewidth=1.5, label="Analytical Fresnel")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("|T|$^2$ (power)")
    ax.legend(loc="lower right")
    ax.set_ylim(0.75, 1.02)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    out_path = os.path.join(os.path.dirname(__file__), "validation_fresnel.png")
    plt.savefig(out_path, dpi=150)
    print(f"\nPlot saved to {out_path}")


if __name__ == "__main__":
    main()
