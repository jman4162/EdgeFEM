#!/usr/bin/env python3
"""openEMS (FDTD) cross-check of the probe-fed patch that EdgeFEM simulates.

Builds the same antenna as StackedPatchDesign / PatchAntennaDesign in
openEMS and reports f_res (|S11| minimum), min |S11| and Z_in. Numbers and
method: docs/validation/patch_2p45ghz.md.

How to run (from the repository root): openEMS runs in the Docker image
`antenna-cad-openems` that antenna-cad builds
(`docker build -t antenna-cad-openems docker/` in the antenna-cad repo; it
has the openEMS Python bindings at /opt/openEMS/venv and works in /sim).
This script copies itself into a run directory, runs itself in the
container with --in-container, and reads the results back. The host side
needs Python 3 with numpy only.

    python3 tools/openems_patch_crosscheck.py --case fr4 \
        --workdir /tmp/oems_fr4 [--edge-res 0.4] [--ground finite]

Output: <workdir>/results.npz (f, s11, zin) and <workdir>/summary.json,
plus a one-line summary on stdout.

Model (units mm; the ground plane spans 0..G in x and y, patch centered):
- Ground: PEC sheet at z = 0, G x G, with a square hole of half-side
  b = probe_hole_radius / (2 * 0.59) around the probe.
- Probe: PEC square post of half-side a = probe_radius / (2 * 0.59) from
  z = 0 to the patch (EdgeFEM's equivalent-radius square post).
- Port: 50 ohm lumped port in the ground-plane hole at z = 0, from the
  probe face x = cx + a to the hole edge x = cx + b, |y - cy| <= a, E
  along +x: the strip on which EdgeFEM puts its gap port.
- Substrate: G x G x h, eps_r with conductivity kappa = 2 pi f0 eps0
  eps_r tan_d, so tan_d is exact at f0 = the case's design frequency.
  (EdgeFEM uses a constant complex eps_r: tan_d at every frequency.)
- Patch: PEC sheet W (x) by L (y) at z = h.
- Air: lambda/4 at the center frequency around the board on every side,
  including below the ground, then 8 PML cells.

--feed vertical replaces the post, hole and gap port with the openEMS
tutorial feed: a z-directed lumped port from the solid ground to the patch.

--ground infinite instead extends the ground and substrate through the
PML, which removes the finite-ground edges; EdgeFEM's domain ends at the
ground edge with an ABC there, so its model lies between the two.

Differences from the EdgeFEM model that remain (see the validation doc):
the hole outside the port strip is open to the air below here and a PMC
(natural) boundary in EdgeFEM; EdgeFEM has no air below or beside the
ground; the loss model above; FDTD staircasing is absent (all faces are
axis-aligned) but the post and gap are resolved by 4 and 3 cells.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time

C0 = 299792458.0
EPS0 = 8.8541878128e-12
DOCKER_IMAGE = "antenna-cad-openems"
SQUARE_EQUIV_RADIUS = 0.59             # as StackedPatchDesign
HOLE_RATIO = math.exp(50.0 / 60.0)     # default probe_hole_radius ratio

# Same cases as tools/patch_convergence.py (dimensions in mm).
CASES = {
    "fr4": dict(L=28.6, W=37.3, h=1.6, eps_r=4.4, tan_d=0.02,
                ground=3 * 37.3, probe=(0.0, -9.0), probe_radius=0.5,
                f0=2.45e9, f_start=1.8e9, f_stop=3.0e9),
    "ro4350b_10g": dict(L=7.6792885, W=9.8200285, h=0.508, eps_r=3.66,
                        tan_d=0.0037, ground=24.809651,
                        probe=(0.0, -7.6792885 / 2 + 2.1307054),
                        probe_radius=0.5, f0=10e9,
                        f_start=8.0e9, f_stop=12.0e9),
}


def _thirds(edge, inside_sign, res):
    """Thirds-rule lines at a metal edge: res/3 inside, 2 res/3 outside."""
    return [edge + inside_sign * res / 3, edge - inside_sign * 2 * res / 3]


def _subdivide(x0, x1, n):
    return [x0 + (x1 - x0) * i / n for i in range(n + 1)]


def run_in_container(params_path):
    import numpy as np
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS

    with open(params_path) as fh:
        prm = json.load(fh)
    c = dict(CASES[prm["case"]])
    if prm.get("probe_radius"):
        c["probe_radius"] = prm["probe_radius"]
    run_dir = os.path.dirname(os.path.abspath(params_path))
    G, L, W, h = c["ground"], c["L"], c["W"], c["h"]
    a = c["probe_radius"] / (2 * SQUARE_EQUIV_RADIUS)
    b = HOLE_RATIO * c["probe_radius"] / (2 * SQUARE_EQUIV_RADIUS)
    cx, cy = G / 2 + c["probe"][0], G / 2 + c["probe"][1]
    x0, x1 = G / 2 - W / 2, G / 2 + W / 2
    y0, y1 = G / 2 - L / 2, G / 2 + L / 2
    f_start, f_stop = c["f_start"], c["f_stop"]
    f_mid = 0.5 * (f_start + f_stop)
    lam_mm = C0 / f_mid * 1e3
    pad = lam_mm / 4
    res = prm["edge_res"] if prm["edge_res"] else h / 4
    # lambda/20 in the substrate at f_stop, lambda/20 in air elsewhere.
    max_res = C0 / f_stop * 1e3 / 20
    n_pml = 8

    fdtd = openEMS(NrTS=prm["max_timesteps"], EndCriteria=prm["end_criteria"])
    fdtd.SetGaussExcite(f_mid, 0.5 * (f_stop - f_start))
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)

    infinite = prm["ground"] == "infinite"
    lo = -pad - n_pml * max_res
    hi = G + pad + n_pml * max_res
    zlo = -pad - n_pml * max_res
    zhi = h + pad + n_pml * max_res
    gx0, gx1 = (lo, hi) if infinite else (0.0, G)

    kappa = 2 * math.pi * c["f0"] * EPS0 * c["eps_r"] * c["tan_d"]
    sub = csx.AddMaterial("substrate", epsilon=c["eps_r"], kappa=kappa)
    sub.AddBox(priority=0, start=[gx0, gx0, 0], stop=[gx1, gx1, h])

    gnd = csx.AddMetal("ground")
    csx.AddMetal("patch").AddBox(priority=10, start=[x0, y0, h],
                                 stop=[x1, y1, h])
    if prm.get("feed", "gap") == "gap":
        for start, stop in (([gx0, gx0], [cx - b, gx1]),
                            ([cx + b, gx0], [gx1, gx1]),
                            ([cx - b, gx0], [cx + b, cy - b]),
                            ([cx - b, cy + b], [cx + b, gx1])):
            gnd.AddBox(priority=10, start=[start[0], start[1], 0],
                       stop=[stop[0], stop[1], 0])
        csx.AddMetal("probe").AddBox(priority=10, start=[cx - a, cy - a, 0],
                                     stop=[cx + a, cy + a, h])
        port = fdtd.AddLumpedPort(1, 50.0, [cx + a, cy - a, 0],
                                  [cx + b, cy + a, 0], "x", 1.0, priority=5)
    else:
        # Tutorial-style feed: solid ground, no post, a z-directed lumped
        # port from the ground to the patch on one mesh line.
        gnd.AddBox(priority=10, start=[gx0, gx0, 0], stop=[gx1, gx1, 0])
        port = fdtd.AddLumpedPort(1, 50.0, [cx, cy, 0], [cx, cy, h], "z",
                                  1.0, priority=5)

    feed_x = (_subdivide(cx - b, cx - a, 3) + _subdivide(cx - a, cx + a, 4)
              + _subdivide(cx + a, cx + b, 3))
    feed_y = (_subdivide(cy - b, cy - a, 3) + _subdivide(cy - a, cy + a, 4)
              + _subdivide(cy + a, cy + b, 3))
    xs = [lo, hi, x0, x1, *feed_x, *_thirds(x0, +1, res),
          *_thirds(x1, -1, res)]
    ys = [lo, hi, y0, y1, *feed_y, *_thirds(y0, +1, res),
          *_thirds(y1, -1, res)]
    if not infinite:
        xs += [0.0, G]
        ys += [0.0, G]
    zs = [zlo, zhi, *_subdivide(0.0, h, prm["n_sub"])]
    for axis, lines in (("x", xs), ("y", ys), ("z", zs)):
        merged = []
        for v in sorted(lines):
            if not merged or v - merged[-1] > 1e-3:
                merged.append(v)
        grid.AddLine(axis, merged)
    grid.SmoothMeshLines("all", max_res, 1.4)
    ncell = [len(grid.GetLines(ax)) for ax in "xyz"]

    t0 = time.time()
    fdtd.Run(os.path.join(run_dir, "fdtd"), cleanup=True)
    t_run = time.time() - t0

    f = np.linspace(f_start, f_stop, prm["n_freq"])
    port.CalcPort(os.path.join(run_dir, "fdtd"), f)
    s11 = port.uf_ref / port.uf_inc
    zin = port.uf_tot / port.if_tot
    np.savez(os.path.join(run_dir, "results.npz"), f=f, s11=s11, zin=zin)
    summary = {"cells": ncell, "t_run_s": t_run, "edge_res_mm": res,
               "max_res_mm": max_res, "kappa": kappa}
    with open(os.path.join(run_dir, "container.json"), "w") as fh:
        json.dump(summary, fh)
    return 0


def _analyse(npz_path, case):
    import numpy as np

    d = np.load(npz_path)
    f, s11, zin = d["f"], d["s11"], d["zin"]
    db = 20 * np.log10(np.abs(s11) + 1e-30)
    # Ignore the outer 5% of the band, where the Gaussian has little energy.
    guard = max(1, len(f) // 20)
    i = guard + int(np.argmin(db[guard:-guard]))
    j = guard + int(np.argmax(zin.real[guard:-guard]))
    out = {
        "f_res_hz": float(f[i]), "min_s11_db": float(db[i]),
        "z_at_res": [float(zin[i].real), float(zin[i].imag)],
        "f_max_re_z_hz": float(f[j]), "max_re_z": float(zin[j].real),
    }
    for ff in (case["f_start"], case["f_stop"]):
        k = int(np.argmin(np.abs(f - ff)))
        out[f"re_z_at_{ff / 1e9:g}ghz"] = float(zin[k].real)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", choices=sorted(CASES), default="fr4")
    ap.add_argument("--workdir", required=False)
    ap.add_argument("--edge-res", type=float, default=None,
                    help="mesh step at the patch edges, mm (default h/4)")
    ap.add_argument("--n-sub", type=int, default=4,
                    help="cells across the substrate thickness")
    ap.add_argument("--ground", choices=["finite", "infinite"],
                    default="finite")
    ap.add_argument("--probe-radius", type=float, default=None,
                    help="probe radius, mm (default: the case's, as EdgeFEM)")
    ap.add_argument("--feed", choices=["gap", "vertical"], default="gap",
                    help="gap: EdgeFEM's post and ground-plane gap port "
                         "(default); vertical: z-directed lumped port from "
                         "ground to patch, no post or hole")
    ap.add_argument("--n-freq", type=int, default=1201)
    ap.add_argument("--end-criteria", type=float, default=1e-5)
    ap.add_argument("--max-timesteps", type=int, default=400000)
    ap.add_argument("--in-container", metavar="PARAMS_JSON")
    args = ap.parse_args()

    if args.in_container:
        return run_in_container(args.in_container)
    if not args.workdir:
        ap.error("--workdir is required")

    run_dir = os.path.abspath(args.workdir)
    os.makedirs(run_dir, exist_ok=True)
    params = {"case": args.case, "edge_res": args.edge_res,
              "n_sub": args.n_sub, "ground": args.ground,
              "probe_radius": args.probe_radius, "feed": args.feed,
              "n_freq": args.n_freq, "end_criteria": args.end_criteria,
              "max_timesteps": args.max_timesteps}
    with open(os.path.join(run_dir, "params.json"), "w") as fh:
        json.dump(params, fh, indent=2)
    script = os.path.join(run_dir, "openems_patch_crosscheck.py")
    shutil.copyfile(os.path.abspath(__file__), script)
    cmd = ["docker", "run", "--rm", "-v", f"{run_dir}:/sim", DOCKER_IMAGE,
           "/sim/openems_patch_crosscheck.py", "--in-container",
           "/sim/params.json"]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    with open(os.path.join(run_dir, "run.log"), "w") as fh:
        fh.write(proc.stdout + "\n" + proc.stderr)
    if proc.returncode != 0:
        print(proc.stderr[-3000:], file=sys.stderr)
        return proc.returncode

    case = CASES[args.case]
    summary = dict(params)
    with open(os.path.join(run_dir, "container.json")) as fh:
        summary.update(json.load(fh))
    summary.update(_analyse(os.path.join(run_dir, "results.npz"), case))
    summary["wall_s"] = wall
    with open(os.path.join(run_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"openEMS {args.case} ground={args.ground} feed={args.feed} "
          f"probe_radius={args.probe_radius or case['probe_radius']} mm "
          f"edge_res={summary['edge_res_mm']:.3g} mm cells={summary['cells']}"
          f" f_res={summary['f_res_hz'] / 1e9:.4f} GHz"
          f" min|S11|={summary['min_s11_db']:.2f} dB"
          f" Z={summary['z_at_res'][0]:.1f}{summary['z_at_res'][1]:+.1f}j"
          f" f(maxReZ)={summary['f_max_re_z_hz'] / 1e9:.4f} GHz"
          f" maxReZ={summary['max_re_z']:.1f} wall={wall:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
