#!/usr/bin/env python3
"""Full-wave FEM sweep of a 2.45 GHz probe-fed patch antenna.

Uses PatchAntennaDesign's FEM backend (lumped port + normalized weights).
Writes S11 vs frequency to patch_fullwave_sparams.csv. Requires gmsh.
"""
from pathlib import Path

import numpy as np

from edgefem.designs import PatchAntennaDesign


def run():
    out = Path(__file__).resolve().parent / "patch_fullwave_sparams.csv"

    # 2.45 GHz patch on FR-4 (dimensions from the analytical design script,
    # examples/patch_antenna_design.py)
    patch = PatchAntennaDesign(
        patch_length=28.6e-3,
        patch_width=37.3e-3,
        substrate_height=1.6e-3,
        substrate_eps_r=4.4,
        substrate_tan_d=0.02,
    )
    patch.set_probe_feed(y_offset=-9.0e-3)
    patch.generate_mesh(density=12, design_freq=2.45e9)

    freqs = np.linspace(2.2e9, 2.7e9, 11)
    with open(out, "w", encoding="utf-8") as f:
        f.write("freq_hz,re_s11,im_s11,abs_s11\n")
        for freq in freqs:
            s11, _ = patch.simulate(freq)
            f.write(f"{freq},{s11.real},{s11.imag},{abs(s11)}\n")
            print(f"{freq/1e9:.3f} GHz  |S11| = {abs(s11):.4f}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    run()
