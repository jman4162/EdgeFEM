#!/usr/bin/env python3
"""
Example 03: Phased Array Unit Cell with Floquet Ports

This example demonstrates how to analyze a phased array unit cell using
periodic boundary conditions with Floquet phase shift. This is the standard
approach for computing embedded element patterns and scan impedance.

Key concepts:
- Periodic boundary conditions pair master/slave edges with phase shift
- Floquet phase shift depends on scan angle and element spacing
- Active S-parameters vary with scan angle

Note: This example uses a simple waveguide as a placeholder for a unit cell.
A real phased array simulation would use a patch antenna or slot geometry.

Requirements:
    - EdgeFEM built with EDGEFEM_PYTHON=ON
    - Mesh file with periodic surfaces

For actual array simulations, the mesh should have:
    - Periodic surfaces in X and/or Y directions
    - A Floquet port (or wave port) for excitation
    - PEC ground plane (if applicable)
"""

import numpy as np
from pathlib import Path

try:
    import pyedgefem as em
except ImportError:
    print("Error: pyedgefem not found. Build EdgeFEM with -DEDGEFEM_PYTHON=ON")
    raise


def compute_floquet_phase(dx, dy, theta_deg, phi_deg, freq):
    """
    Compute Floquet phase shift for given scan angle.

    Parameters:
        dx: Element spacing in x (meters)
        dy: Element spacing in y (meters)
        theta_deg: Scan angle from broadside (degrees)
        phi_deg: Azimuth angle (degrees)
        freq: Frequency (Hz)

    Returns:
        (phase_x, phase_y): Phase shifts in radians
    """
    c0 = 299792458
    k0 = 2 * np.pi * freq / c0

    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)

    # Transverse wave vector components
    kx = k0 * np.sin(theta) * np.cos(phi)
    ky = k0 * np.sin(theta) * np.sin(phi)

    phase_x = kx * dx
    phase_y = ky * dy

    return phase_x, phase_y


def main():
    """Demonstrate Floquet/periodic BC concepts."""

    print("=" * 60)
    print("EdgeFEM Example 03: Phased Array Unit Cell Concepts")
    print("=" * 60)

    # ========================================================================
    # Part 1: Floquet Phase Calculation
    # ========================================================================
    print("\n1. Floquet Phase Shift Calculation")
    print("   " + "-" * 50)

    # Example: half-wavelength spacing at 10 GHz
    freq = 10e9
    c0 = 299792458
    wavelength = c0 / freq
    dx = wavelength / 2  # ~15 mm at 10 GHz
    dy = wavelength / 2

    print(f"   Frequency: {freq/1e9:.1f} GHz")
    print(f"   Wavelength: {wavelength*1000:.1f} mm")
    print(f"   Element spacing: dx = dy = {dx*1000:.1f} mm (lambda/2)")

    # Compute phase shift for various scan angles
    print(f"\n   Scan angle vs Floquet phase:")
    print(f"   {'Theta(deg)':<12} {'Phi(deg)':<10} {'Phase_x(deg)':<14} {'Phase_y(deg)':<14}")
    print("   " + "-" * 50)

    scan_angles = [
        (0, 0),     # Broadside
        (15, 0),    # 15 deg in E-plane
        (30, 0),    # 30 deg in E-plane
        (45, 0),    # 45 deg in E-plane
        (30, 45),   # D-plane
        (30, 90),   # H-plane
    ]

    for theta, phi in scan_angles:
        px, py = compute_floquet_phase(dx, dy, theta, phi, freq)
        print(f"   {theta:<12} {phi:<10} {np.degrees(px):<14.1f} {np.degrees(py):<14.1f}")

    # ========================================================================
    # Part 2: Active Impedance Concepts
    # ========================================================================
    print("\n2. Active Impedance Concepts")
    print("   " + "-" * 50)

    print("""
   For a phased array, the active impedance at element n is:

       Z_active,n = Z0 * (1 + Gamma_active,n) / (1 - Gamma_active,n)

   where the active reflection coefficient is:

       Gamma_active,n = sum_m S_nm * a_m / a_n

   For uniform excitation with progressive phase:
       a_m = exp(j * m * psi)

   where psi is the inter-element phase shift.

   EdgeFEM provides these functions:
   - compute_active_impedance(coupling, excitation)
   - active_impedance_scan(coupling, positions, theta, phi, k0)
   - active_vswr(Gamma_active)
   """)

    # ========================================================================
    # Part 3: Using Periodic BCs in EdgeFEM
    # ========================================================================
    print("\n3. Using Periodic BCs in EdgeFEM")
    print("   " + "-" * 50)

    print("""
   For actual unit cell simulations, the workflow is:

   1. Create mesh with periodic surfaces (master/slave pairs)

   2. Build periodic boundary conditions:
      pbc = em.build_periodic_pairs(mesh, master_tag, slave_tag,
                                     period_vector, tolerance)

   3. Set Floquet phase shift:
      k_transverse = [kx, ky, 0]  # From scan angle
      em.set_floquet_phase(pbc, k_transverse)

   4. Assemble and solve with periodic assembly:
      asmbl = em.assemble_maxwell_periodic(mesh, params, bc, pbc, ports, 0)
      result = em.solve_linear(asmbl.A, asmbl.b, options)

   5. Extract S-parameters with periodic BC:
      S = em.calculate_sparams_periodic(mesh, params, bc, pbc, ports)

   Example code:

   ```python
   # Build periodic pairs
   period_x = np.array([dx, 0, 0])
   pbc = em.build_periodic_pairs(mesh, x_minus_tag, x_plus_tag,
                                  period_x, 1e-9)

   # Set phase for 30 degree scan in E-plane
   k0 = 2 * np.pi * freq / c0
   theta = np.radians(30)
   k_trans = np.array([k0 * np.sin(theta), 0, 0])
   em.set_floquet_phase(pbc, k_trans)

   # Calculate S-parameters
   S = em.calculate_sparams_periodic(mesh, params, bc, pbc, ports)
   ```
   """)

    # ========================================================================
    # Part 4: Scan Loss Estimation
    # ========================================================================
    print("\n4. Scan Loss vs Scan Angle")
    print("   " + "-" * 50)

    # Theoretical cos(theta) pattern
    print("   Ideal cos(theta) scan loss:")
    print(f"   {'Theta(deg)':<12} {'Scan Loss (dB)':<15}")
    print("   " + "-" * 30)

    for theta in [0, 15, 30, 45, 60]:
        if theta == 0:
            scan_loss = 0
        else:
            scan_loss = -20 * np.log10(np.cos(np.radians(theta)))
        print(f"   {theta:<12} {scan_loss:<15.2f}")

    print("\n   Note: Actual scan loss includes:")
    print("   - Element pattern (gain reduction at angle)")
    print("   - Scan impedance mismatch")
    print("   - Grating lobe onset (if spacing > lambda)")
    print("   - Scan blindness (surface wave coupling)")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("""
This example demonstrated:
1. Floquet phase shift calculation for phased arrays
2. Active impedance concepts for embedded element analysis
3. EdgeFEM periodic BC workflow for unit cell simulations
4. Scan loss estimation for phased array design

For a complete unit cell simulation, you would need:
- A mesh with the actual antenna geometry (patch, slot, etc.)
- Periodic surfaces matching the array lattice
- Floquet or wave ports for excitation

EdgeFEM provides all the necessary functions for full-wave
unit cell analysis of phased array antennas.
""")

    print("Demo complete.")


if __name__ == "__main__":
    main()
