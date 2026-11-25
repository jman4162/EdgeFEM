#!/usr/bin/env python3
"""Periodic unit cell demo: reflection coefficient of a bare grounded cell.

EXPERIMENTAL, NOT VALIDATED. UnitCellDesign has no validated Floquet port:
its port uses the 1.0.0 modal terms (on this bare grounded slab |R| = 1 but
arg(R) is 30-106 deg off the analytic value), the periodic BC covers the x
faces only, and patches are not meshed. The results below need
experimental=True and must not be used for design.
"""

import numpy as np

from edgefem.designs import UnitCellDesign


def main():
    # 5 mm square cell, thin grounded substrate (no patch: not meshed yet)
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=0.5e-3,
        substrate_eps_r=3.5,
        experimental=True,
    )

    print("Generating mesh (coarse)...")
    cell.generate_mesh(density=8)

    print("\nReflection vs frequency at normal incidence:")
    print("  Freq (GHz) |   |R|    | phase(R) (deg)")
    print("  " + "-" * 40)
    for freq in np.linspace(8e9, 12e9, 5):
        R, T = cell.reflection_transmission(freq, theta=0.0, phi=0.0)
        print(f"  {freq/1e9:9.2f} | {abs(R):7.4f} | {np.degrees(np.angle(R)):+9.1f}")

    Zs = cell.surface_impedance(10e9)
    print(f"\nSurface impedance at 10 GHz: {Zs:.1f} ohm")


if __name__ == "__main__":
    main()
