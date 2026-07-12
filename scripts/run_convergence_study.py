#!/usr/bin/env python3
"""WR-90 mesh convergence study.

Generates WR-90 waveguide meshes at several densities, runs the full-wave
eigenmode S-parameter solve at 10 GHz via RectWaveguideDesign, compares
against the analytical TE10 result (|S21| = 1, phase = -beta*L), and fits
the observed convergence order from the |S21| error.

Regenerates the convergence table in docs/validation.md:

    PYTHONPATH=build/python:python python3 scripts/run_convergence_study.py

Requires gmsh and an EdgeFEM build with EDGEFEM_PYTHON=ON.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from edgefem.designs import RectWaveguideDesign

# WR-90 at 10 GHz
A = 22.86e-3
B = 10.16e-3
LENGTH = 50e-3
FREQ = 10e9

C0 = 299792458.0


def analytical_phase_deg(freq: float) -> float:
    k0 = 2 * np.pi * freq / C0
    kc = np.pi / A
    beta = np.sqrt(k0**2 - kc**2)
    return np.degrees(-beta * LENGTH)


def wrap_deg(angle: float) -> float:
    """Wrap an angle difference to [-180, 180) degrees."""
    return (angle + 180.0) % 360.0 - 180.0


def run_point(density: float, out_dir: Path) -> dict:
    wg = RectWaveguideDesign(a=A, b=B, length=LENGTH)
    mesh_path = str(out_dir / f"wr90_d{density:g}.msh")
    wg.generate_mesh(density=density, output_path=mesh_path)

    t0 = time.perf_counter()
    S = wg.sparams_at_freq(FREQ)
    runtime = time.perf_counter() - t0

    s21 = S[1, 0]
    phase_err = abs(wrap_deg(np.degrees(np.angle(s21)) - analytical_phase_deg(FREQ)))
    return {
        "density": density,
        "n_tets": len(wg.mesh.tets),
        "n_edges": len(wg.mesh.edges),
        "abs_s11": float(abs(S[0, 0])),
        "abs_s21": float(abs(s21)),
        "s21_error": float(abs(abs(s21) - 1.0)),
        "phase_error_deg": float(phase_err),
        "runtime_s": runtime,
    }


def observed_order(results: list) -> float:
    """Least-squares slope of log(|S21| error) vs log(h), h ~ density^-1."""
    h = np.array([1.0 / r["density"] for r in results])
    err = np.array([max(r["s21_error"], 1e-12) for r in results])
    slope, _ = np.polyfit(np.log(h), np.log(err), 1)
    return float(slope)


def main():
    parser = argparse.ArgumentParser(description="WR-90 mesh convergence study")
    parser.add_argument("--output-dir", default="convergence_data")
    parser.add_argument("--densities", type=float, nargs="+",
                        default=[5, 8, 10, 14, 18])
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    print(f"WR-90 convergence study at {FREQ/1e9:.0f} GHz "
          f"(analytical phase {analytical_phase_deg(FREQ):.1f} deg)\n")

    results = []
    for density in args.densities:
        print(f"density = {density:g} elements/wavelength ...", flush=True)
        r = run_point(density, out_dir)
        results.append(r)
        print(f"  tets={r['n_tets']}, |S21|={r['abs_s21']:.4f}, "
              f"|S11|={r['abs_s11']:.4f}, phase err={r['phase_error_deg']:.1f} deg, "
              f"{r['runtime_s']:.1f} s")

    order = observed_order(results)

    with open(out_dir / "convergence_results.json", "w") as f:
        json.dump({"freq_hz": FREQ, "observed_order": order,
                   "results": results}, f, indent=2)

    print("\n| Elements/λ | Tets | |S21| | |S21| error | |S11| | Phase err | Runtime |")
    print("|-----------|------|-------|-------------|-------|-----------|---------|")
    for r in results:
        print(f"| {r['density']:g} | {r['n_tets']:,} | {r['abs_s21']:.4f} "
              f"| {r['s21_error']*100:.2f}% | {r['abs_s11']:.4f} "
              f"| {r['phase_error_deg']:.1f}° | {r['runtime_s']:.1f}s |")
    print(f"\nObserved |S21|-error convergence order: {order:.2f} "
          f"(w.r.t. h ~ 1/density)")
    print(f"Results written to {out_dir}/convergence_results.json")


if __name__ == "__main__":
    main()
