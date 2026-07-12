#!/usr/bin/env python3
"""Periodic unit cell demo: reflection coefficient of a patch-loaded cell.

Runs a real periodic FEM solve via UnitCellDesign.reflection_transmission().

Scope (see README "Current Limitations"): the solver applies the Bloch phase
along the x lattice axis only, the port is the cell's dominant waveguide mode
(not a Floquet harmonic), and the periodic constraint elimination is dense,
so keep meshes small. Results are most meaningful at normal incidence.
"""

import numpy as np

from edgefem.designs import UnitCellDesign


def main():
    # 5 mm square cell, patch on a thin grounded substrate
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=0.5e-3,
        substrate_eps_r=3.5,
    )
    cell.add_patch(width=4e-3, length=4e-3)

    print("Generating mesh (coarse; periodic elimination is dense)...")
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
