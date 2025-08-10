#!/usr/bin/env python3
"""Analytic rectangular microstrip patch antenna design with simple sweeps.

This script implements closed-form design equations for a rectangular
microstrip patch antenna on an FR-4 substrate. It produces quick parametric
sweeps for the inset depth and patch length and exports S-parameter data
using a simple RLC resonance model.

The implementation is intentionally lightweight and does not call a full EM
solver. Results are approximate but follow the design targets in
``patch_antenna_design_prompt.md``.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from math import acos, cos, pi, sqrt
from typing import Iterable, Tuple

import numpy as np

C0 = 299_792_458.0  # speed of light (m/s)


@dataclass
class Substrate:
    er: float
    h: float  # thickness (m)


@dataclass
class PatchDesign:
    freq: float  # design frequency (Hz)
    substrate: Substrate

    def patch_dimensions(self) -> Tuple[float, float, float, float]:
        """Return (W, L, eps_eff, delta_L)."""
        er, h = self.substrate.er, self.substrate.h
        f = self.freq
        W = C0 / (2 * f) * sqrt(2 / (er + 1))
        eps_eff = (er + 1) / 2 + (er - 1) / (2 * sqrt(1 + 12 * h / W))
        F1 = (eps_eff + 0.3) * (W / h + 0.264)
        F2 = (eps_eff - 0.258) * (W / h + 0.8)
        delta_L = 0.412 * h * F1 / F2
        L = C0 / (2 * f * sqrt(eps_eff)) - 2 * delta_L
        return W, L, eps_eff, delta_L

    @staticmethod
    def edge_resistance(W: float, L: float, er: float, eps_eff: float) -> float:
        return 90 * (eps_eff ** 2 / (er - 1)) * (L / W) ** 2

    @staticmethod
    def inset_resistance(R_edge: float, L: float, inset: float) -> float:
        # Feeding deeper into the patch decreases the input resistance.
        # The edge resistance is scaled by cos^2 of the inset position.
        return R_edge * (cos(pi * inset / L) ** 2)

    @staticmethod
    def resonance_freq(L: float, eps_eff: float, delta_L: float) -> float:
        return C0 / (2 * (L + 2 * delta_L) * sqrt(eps_eff))

    @staticmethod
    def input_impedance(f: float, f0: float, R: float, Q: float) -> complex:
        X = R * 2 * Q * (f / f0 - f0 / f)
        return complex(R, X)


def s11_from_z(z: complex) -> complex:
    return (z - 50) / (z + 50)


def sweep_inset(design: PatchDesign, y0_range: Iterable[float], W: float, L: float,
                R_edge: float) -> Tuple[np.ndarray, np.ndarray]:
    s11 = []
    for y0 in y0_range:
        Rin = design.inset_resistance(R_edge, L, y0)
        gamma = s11_from_z(complex(Rin, 0))
        s11.append(20 * np.log10(abs(gamma)))
    return np.array(y0_range), np.array(s11)


def sweep_length(design: PatchDesign, L_vals: Iterable[float], inset: float,
                 W: float, eps_eff: float, delta_L: float, Q: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    fres, s11 = [], []
    for L in L_vals:
        f0 = design.resonance_freq(L, eps_eff, delta_L)
        R_edge = design.edge_resistance(W, L, design.substrate.er, eps_eff)
        Rin = design.inset_resistance(R_edge, L, inset)
        Zin = design.input_impedance(design.freq, f0, Rin, Q)
        gamma = s11_from_z(Zin)
        fres.append(f0)
        s11.append(20 * np.log10(abs(gamma)))
    return np.array(L_vals), np.array(fres), np.array(s11)


def export_sparameters(filename: str, design: PatchDesign, freqs: np.ndarray,
                       L: float, inset: float, W: float,
                       eps_eff: float, delta_L: float, Q: float) -> None:
    R_edge = design.edge_resistance(W, L, design.substrate.er, eps_eff)
    Rin = design.inset_resistance(R_edge, L, inset)
    f0 = design.resonance_freq(L, eps_eff, delta_L)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("# Hz S RI R 50\n")
        for fr in freqs:
            Zin = design.input_impedance(fr, f0, Rin, Q)
            gamma = s11_from_z(Zin)
            f.write(f"{fr} {gamma.real} {gamma.imag}\n")


def export_csv(filename: str, design: PatchDesign, freqs: np.ndarray,
               L: float, inset: float, W: float,
               eps_eff: float, delta_L: float, Q: float) -> None:
    R_edge = design.edge_resistance(W, L, design.substrate.er, eps_eff)
    Rin = design.inset_resistance(R_edge, L, inset)
    f0 = design.resonance_freq(L, eps_eff, delta_L)

    with open(filename, "w", encoding="utf-8") as f:
        f.write("freq_hz,s11_db,rin,xin\n")
        for fr in freqs:
            Zin = design.input_impedance(fr, f0, Rin, Q)
            gamma = s11_from_z(Zin)
            f.write(f"{fr},{20*np.log10(abs(gamma))},{Zin.real},{Zin.imag}\n")


def export_geo(filename: str, W: float, L: float, inset: float, feed_w: float,
               board_x: float = 0.08, board_y: float = 0.07,
               h: float = 1.6e-3, air: float = 0.0306) -> None:
    """Write a minimal parametric Gmsh geometry file."""
    with open(filename, "w", encoding="utf-8") as g:
        g.write('SetFactory("OpenCASCADE");\n')
        g.write('// Units: meters\n')
        g.write(f'h={h}; W={W}; L={L}; inset={inset}; feedW={feed_w};\n')
        g.write(f'board_x={board_x}; board_y={board_y}; air={air};\n')
        g.write('Rectangle(1) = {-board_x/2, -board_y/2, 0, board_x, board_y};\n')
        g.write('Rectangle(2) = {-W/2, -L/2, h, W, L};\n')
        g.write('Rectangle(3) = {-feedW/2, -board_y/2, h, feedW, board_y/2 - L/2 + inset};\n')
        g.write('Rectangle(4) = {-feedW/2, -L/2, h, feedW, inset};\n')
        g.write('tmp1[] = BooleanDifference{ Surface{2}; Delete; }{ Surface{4}; Delete; };\n')
        g.write('tmp2[] = BooleanUnion{ Surface{tmp1[0]}; Delete; }{ Surface{3}; Delete; };\n')
        g.write('Rectangle(5) = {-(board_x/2+air), -(board_y/2+air), -air, board_x+2*air, board_y+2*air};\n')
        g.write('Extrude {0,0,h} { Surface{1}; }\n')
        g.write('Extrude {0,0,h+air} { Surface{5}; }\n')


def main() -> None:
    sub = Substrate(er=4.4, h=1.6e-3)
    design = PatchDesign(freq=2.45e9, substrate=sub)
    W, L, eps_eff, delta_L = design.patch_dimensions()
    R_edge = design.edge_resistance(W, L, sub.er, eps_eff)

    # Inset depth sweep 4..11 mm
    y0_vals_m = np.linspace(4e-3, 11e-3, 8)
    y0_vals, s11_y0 = sweep_inset(design, y0_vals_m, W, L, R_edge)
    best_idx = np.argmin(s11_y0)
    best_y0 = y0_vals_m[best_idx]

    # Length sweep around nominal L
    L_vals = np.linspace(L - 0.6e-3, L + 0.6e-3, 7)
    Q = design.freq / 80e6  # target ~80 MHz bandwidth
    L_sweep, fres, s11_L = sweep_length(design, L_vals, best_y0, W, eps_eff, delta_L, Q)
    best_L = L_vals[np.argmin(s11_L)]

    # Frequency sweep for final design
    freqs = np.linspace(2.2e9, 2.7e9, 101)
    export_sparameters("examples/patch_2p45_FR4.s1p", design, freqs, best_L, best_y0, W, eps_eff, delta_L, Q)
    export_csv("examples/patch_sparams.csv", design, freqs, best_L, best_y0, W, eps_eff, delta_L, Q)
    export_geo("examples/patch_antenna.geo", W, best_L, best_y0, 3.08e-3)

    summary = {
        "freq_hz": design.freq,
        "substrate": {"er": sub.er, "h_m": sub.h},
        "patch_width_m": W,
        "patch_length_m": best_L,
        "inset_depth_m": best_y0,
        "feed_width_m": 3.08e-3,
        "quality_factor": Q,
    }
    with open("examples/patch_design.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    # Console report
    print("Inset depth sweep (mm vs S11 dB):")
    for mm_val, s in zip(y0_vals * 1e3, s11_y0):
        print(f"  {mm_val:4.1f} mm : {s:6.2f} dB")
    print(f"Best inset ≈ {best_y0*1e3:.2f} mm")

    print("Length sweep (mm, f_res GHz, S11 dB):")
    for Lm, fr, s in zip(L_sweep * 1e3, fres, s11_L):
        print(f"  {Lm:5.2f} mm : {fr/1e9:6.3f} GHz : {s:6.2f} dB")
    print(f"Best length ≈ {best_L*1e3:.2f} mm")

    print(
        "Outputs written: examples/patch_sparams.csv, examples/patch_2p45_FR4.s1p,"
        " examples/patch_antenna.geo, examples/patch_design.json"
    )


if __name__ == "__main__":
    main()
