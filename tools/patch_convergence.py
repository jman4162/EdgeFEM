#!/usr/bin/env python3
"""Mesh convergence study for the probe-fed patch (StackedPatchDesign).

Reproduces the convergence table in docs/validation/patch_2p45ghz.md.
Each (case, density, patch factor) runs in its own subprocess so the peak
resident memory of one run can be read from getrusage; runs go one at a
time (a density-12 FR-4 run peaks near 6.5 GB).

Run from the repository root with a build that has pyedgefem:

    PYTHONPATH=build-fix/python:python python tools/patch_convergence.py \
        --case fr4 --density 8 12 16 --out conv.jsonl

    # near-patch mesh size at density 12: 6, 4 and 3 x h / density
    PYTHONPATH=build-fix/python:python python tools/patch_convergence.py \
        --case fr4 --density 12 --patch-factor 6 4 3 --out conv.jsonl

    # EdgeArray's 10 GHz RO4350B patch
    PYTHONPATH=build-fix/python:python python tools/patch_convergence.py \
        --case ro4350b_10g --density 6 12 --out ro.jsonl

Cases:
  fr4          2.45 GHz FR-4 patch of the known-answer test
               (python/tests/test_patch_physics.py): L 28.6 mm, W 37.3 mm,
               h 1.6 mm, eps_r 4.4, tan d 0.02, ground 111.9 mm, probe
               radius 0.5 mm at y = -9 mm. Air box lambda/2 at 1.8 GHz.
               Sweeps 1.8-3.0 GHz (100 MHz) and 2.25-2.60 GHz (10 MHz).
  ro4350b_10g  EdgeArray's synthesized 10 GHz patch on 0.508 mm RO4350B
               (L 7.679 mm, W 9.820 mm, eps_r 3.66, tan d 0.0037, ground
               24.81 mm), probe radius 0.5 mm at the inset point,
               y = -L/2 + 2.131 mm, as EdgeArray's EdgeFEM worker builds it.
               Air box lambda/2 at 8.5 GHz. Sweeps 8.5-11.5 GHz (100 MHz,
               EdgeArray's live-test sweep) and 9.0-10.5 GHz (25 MHz).

Each run prints one JSON line: mesh size, wall times, peak memory, the
|S11| minimum of the fine sweep (grid point and parabolic refinement in
dB), Z there, and Re(Z) at the ends of the wide sweep.
"""

import argparse
import json
import math
import os
import resource
import subprocess
import sys
import tempfile
import time

C0 = 299792458.0

CASES = {
    "fr4": dict(
        patch=dict(patch_length=28.6e-3, patch_width=37.3e-3,
                   substrate_height=1.6e-3, substrate_eps_r=4.4,
                   substrate_tan_d=0.02),
        ground=None,  # PatchAntennaDesign default, 3 x 37.3 mm
        probe=(0.0, -9.0e-3),
        design_freq=2.45e9,
        f_min=1.8e9,
        wide=(1.8e9, 3.0e9, 13),
        fine=(2.25e9, 2.60e9, 36),
    ),
    "ro4350b_10g": dict(
        patch=dict(patch_length=7.6792885e-3, patch_width=9.8200285e-3,
                   substrate_height=0.508e-3, substrate_eps_r=3.66,
                   substrate_tan_d=0.0037),
        ground=24.809651e-3,
        probe=(0.0, -7.6792885e-3 / 2 + 2.1307054e-3),
        design_freq=10e9,
        f_min=8.5e9,
        wide=(8.5e9, 11.5e9, 31),
        fine=(9.0e9, 10.5e9, 61),
    ),
}


def cavity_model_resonance(length, width, h, eps_r):
    """TM010 resonance, Balanis Sec. 14.2 with Hammerstad's dL (as in
    python/tests/test_patch_physics.py)."""
    eps_reff = (eps_r + 1) / 2 + (eps_r - 1) / 2 / math.sqrt(1 + 12 * h / width)
    dl = (0.412 * h * (eps_reff + 0.3) * (width / h + 0.264)
          / ((eps_reff - 0.258) * (width / h + 0.8)))
    return C0 / (2 * (length + 2 * dl) * math.sqrt(eps_reff))


def _parabolic_min(freqs, mag_db, i):
    """Vertex of the parabola through the three points around index i."""
    if i == 0 or i == len(freqs) - 1:
        return float(freqs[i]), float(mag_db[i])
    y0, y1, y2 = mag_db[i - 1], mag_db[i], mag_db[i + 1]
    denom = y0 - 2 * y1 + y2
    if denom <= 0:
        return float(freqs[i]), float(y1)
    step = freqs[i + 1] - freqs[i]
    d = 0.5 * (y0 - y2) / denom
    return float(freqs[i] + d * step), float(y1 - 0.25 * (y0 - y2) * d)


def run_one(case_name, density, patch_factor, air_box_height=None):
    import numpy as np
    from edgefem.designs import PatchAntennaDesign, StackedPatchDesign

    case = CASES[case_name]
    if patch_factor is not None:
        StackedPatchDesign._PATCH_SIZE_X_DENSITY = patch_factor
    kwargs = dict(case["patch"])
    if case["ground"] is not None:
        kwargs["ground_plane_size"] = case["ground"]
    if air_box_height is not None:
        kwargs["air_box_height"] = air_box_height
    patch = PatchAntennaDesign(**kwargs)
    patch.set_probe_feed(x_offset=case["probe"][0], y_offset=case["probe"][1])

    out = {"case": case_name, "density": density,
           "patch_factor": StackedPatchDesign._PATCH_SIZE_X_DENSITY,
           "air_box_height_m": air_box_height}
    with tempfile.TemporaryDirectory() as tmp:
        t0 = time.perf_counter()
        patch.generate_mesh(density=density, design_freq=case["design_freq"],
                            f_min=case["f_min"],
                            output_path=os.path.join(tmp, "patch.msh"))
        out["t_mesh_s"] = time.perf_counter() - t0
    mesh = patch.mesh
    out["tets"] = len(mesh.tet_phys())
    out["edges"] = int(mesh.num_edges())

    t0 = time.perf_counter()
    fw, zw, sw = patch.frequency_sweep(*case["wide"])
    out["t_wide_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    ff, zf, sf = patch.frequency_sweep(*case["fine"])
    out["t_fine_s"] = time.perf_counter() - t0

    # A failed solve gives S11 = NaN (and Z = 1e6 from the design class).
    out["failed_points"] = int(np.isnan(np.asarray(sw)).sum()
                               + np.isnan(np.asarray(sf)).sum())
    if out["failed_points"]:
        return out
    mag_db = 20 * np.log10(np.abs(np.asarray(sf)) + 1e-30)
    i = int(np.argmin(mag_db))
    f_par, db_par = _parabolic_min(np.asarray(ff), mag_db, i)
    zf = np.asarray(zf)
    p = case["patch"]
    f_cav = cavity_model_resonance(p["patch_length"], p["patch_width"],
                                   p["substrate_height"], p["substrate_eps_r"])
    wide_db = 20 * np.log10(np.abs(np.asarray(sw)) + 1e-30)
    iw = int(np.argmin(wide_db))
    phase = np.degrees(np.unwrap(np.angle(np.asarray(sw))))
    out.update({
        "f_cavity_hz": f_cav,
        "f_res_grid_hz": float(ff[i]),
        "f_res_interp_hz": f_par,
        "f_res_vs_cavity_pct": 100 * (f_par - f_cav) / f_cav,
        "min_s11_db_grid": float(mag_db[i]),
        "min_s11_db_interp": db_par,
        "min_at_fine_edge": i in (0, len(ff) - 1),
        "z_at_res": [float(zf[i].real), float(zf[i].imag)],
        "max_re_z_fine": float(zf.real.max()),
        "f_max_re_z_hz": float(ff[int(np.argmax(zf.real))]),
        "wide": {
            "freq_hz": [float(f) for f in fw],
            "z_re": [float(z.real) for z in np.asarray(zw)],
            "z_im": [float(z.imag) for z in np.asarray(zw)],
            "s11_db": [float(d) for d in wide_db],
            "s11_phase_deg": [float(d) for d in phase],
        },
        "wide_min_s11_db": float(wide_db[iw]),
        "wide_f_min_s11_hz": float(fw[iw]),
        "wide_argmin_interior": 0 < iw < len(fw) - 1,
        "wide_phase_spread_deg": float(phase.max() - phase.min()),
        "fine": {
            "freq_hz": [float(f) for f in ff],
            "z_re": [float(z.real) for z in zf],
            "z_im": [float(z.imag) for z in zf],
            "s11_db": [float(d) for d in mag_db],
        },
    })
    # ru_maxrss is bytes on macOS, kilobytes on Linux.
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    out["peak_rss_gb"] = rss / (1e9 if sys.platform == "darwin" else 1e6)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", choices=sorted(CASES), default="fr4")
    ap.add_argument("--density", type=float, nargs="+", default=[8, 12, 16])
    ap.add_argument("--patch-factor", type=float, nargs="+", default=[None],
                    help="override StackedPatchDesign._PATCH_SIZE_X_DENSITY")
    ap.add_argument("--air-box-height", type=float,
                    help="air box height above the patch (m); default "
                         "lambda/2 at the case's f_min")
    ap.add_argument("--out", help="append JSON lines here")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker:
        res = run_one(args.case, args.density[0], args.patch_factor[0],
                      args.air_box_height)
        print(json.dumps(res))
        return 0

    for density in args.density:
        for factor in args.patch_factor:
            cmd = [sys.executable, os.path.abspath(__file__), "--worker",
                   "--case", args.case, "--density", str(density)]
            if factor is not None:
                cmd += ["--patch-factor", str(factor)]
            if args.air_box_height is not None:
                cmd += ["--air-box-height", str(args.air_box_height)]
            t0 = time.perf_counter()
            proc = subprocess.run(cmd, capture_output=True, text=True)
            wall = time.perf_counter() - t0
            if proc.returncode != 0:
                print(f"density {density} factor {factor}: failed "
                      f"(exit {proc.returncode})\n{proc.stderr[-2000:]}",
                      file=sys.stderr)
                continue
            line = proc.stdout.strip().splitlines()[-1]
            res = json.loads(line)
            res["wall_s"] = wall
            if args.out:
                with open(args.out, "a") as fh:
                    fh.write(json.dumps(res) + "\n")
            if res["failed_points"]:
                print(f"{res['case']} d={density:g} "
                      f"factor={res['patch_factor']:g} edges={res['edges']}: "
                      f"{res['failed_points']} solves failed (no field)",
                      flush=True)
                continue
            w = res["wide"]
            print(f"{res['case']} d={density:g} factor={res['patch_factor']:g}"
                  f" tets={res['tets']} edges={res['edges']}"
                  f" f_res={res['f_res_interp_hz'] / 1e9:.4f} GHz"
                  f" ({res['f_res_vs_cavity_pct']:+.2f}% vs cavity"
                  f" {res['f_cavity_hz'] / 1e9:.4f})"
                  f" min|S11|={res['min_s11_db_interp']:.2f} dB"
                  f" Z={res['z_at_res'][0]:.1f}{res['z_at_res'][1]:+.1f}j"
                  f" ReZ[{w['freq_hz'][0] / 1e9:g}]={w['z_re'][0]:.2f}"
                  f" ReZ[{w['freq_hz'][-1] / 1e9:g}]={w['z_re'][-1]:.2f}"
                  f" wall={wall:.0f}s rss={res['peak_rss_gb']:.1f}GB",
                  flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
