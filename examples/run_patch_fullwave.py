#!/usr/bin/env python3
import subprocess
from pathlib import Path
import numpy as np
import pyvectorem as em


def run():
    script_dir = Path(__file__).resolve().parent
    geo = script_dir / "patch_antenna.geo"
    msh = script_dir / "patch_antenna.msh"
    out = script_dir / "patch_fullwave_sparams.csv"
    subprocess.run(["gmsh", str(geo), "-3", "-o", str(msh)], check=True)
    mesh = em.load_gmsh(str(msh))
    bc = em.BC()
    p = em.MaxwellParams()
    freqs = np.linspace(2.2e9, 2.7e9, 11)
    with open(out, "w", encoding="utf-8") as f:
        f.write("freq_hz,residual\n")
        for freq in freqs:
            p.omega = 2 * np.pi * freq
            asm = em.assemble_maxwell(mesh, p, bc)
            res = em.solve_linear(asm.A, asm.b)
            f.write(f"{freq},{abs(res.residual)}\n")


if __name__ == "__main__":
    run()
