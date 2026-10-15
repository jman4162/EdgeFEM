#!/usr/bin/env python3
"""Smoke test for an installed edgefem wheel; runs in a few seconds.

Checks the package version against the wheel metadata, then solves the WR-90
waveguide of examples/rect_waveguide.msh (committed mesh, no gmsh needed) at
10 GHz through `RectWaveguideDesign`. A matched through line must
have |S21| near 1 and |S11| near 0.

    python tools/wheel_smoke.py [--mesh PATH] [--expect-umfpack on|off]

Run it from a directory other than the repository root, or the repository's
python/ directory can shadow the installed package.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--mesh",
        default=str(Path(__file__).resolve().parent.parent / "examples" / "rect_waveguide.msh"),
    )
    ap.add_argument(
        "--expect-umfpack",
        choices=["on", "off"],
        default="off",
        help="expected build: published wheels are built without UMFPACK",
    )
    args = ap.parse_args()

    import edgefem
    import pyedgefem as em

    meta = importlib.metadata.version("edgefem")
    assert edgefem.__version__ == meta, (edgefem.__version__, meta)
    assert em.umfpack_available() == (args.expect_umfpack == "on"), (
        "umfpack_available() =",
        em.umfpack_available(),
    )

    from edgefem.designs.waveguide import RectWaveguideDesign

    wg = RectWaveguideDesign(a=0.02286, b=0.01016, length=0.05)  # WR-90
    wg.load_mesh(args.mesh)
    s = wg.sparams_at_freq(10e9)
    s11, s21 = abs(s[0, 0]), abs(s[1, 0])
    print(f"edgefem {edgefem.__version__}, UMFPACK {em.umfpack_available()}: "
          f"|S11| = {s11:.4f}, |S21| = {s21:.4f}")
    assert s11 < 0.1 and s21 > 0.95, (s11, s21)
    return 0


if __name__ == "__main__":
    sys.exit(main())
