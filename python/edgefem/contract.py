"""Producer-owned artifact contract fixtures.

EdgeFEM owns the ``EdgeFEM_ArrayPackage`` JSON format (written by
``export_array_package_json``) and the embedded-pattern CSV format (written by
``export_patterns_csv``). Consumers — opensatcom's ``edgefem_json_loader`` is
the first — vendor the golden fixture this module writes and pin their loaders
against it, so a format change here breaks a test there instead of silently
producing wrong physics.

The fixture values are canonical hand-pinned numbers, not solver output: the
contract under test is the *format* (keys, shapes, units, conventions), and a
deterministic fixture keeps the test hermetic. Conventions the fixture
demonstrates and the docs guarantee:

- ``element_positions``: metres, always 3-vectors.
- ``frequencies``: Hz, ascending.
- ``S_matrices``: one (N, N) complex matrix per frequency, entries serialized
  as ``[re, im]`` pairs, 50-ohm reference unless ``reference_impedance`` says
  otherwise. These are scattering parameters: a consumer building a coupling
  correction should use the (I + S)^-1 form.
- Embedded patterns (CSV, because the JSON exporter writes only a summary):
  two principal-plane cuts per port, phi = 0 (E-plane) and phi = 90 deg
  (H-plane), columns ``port_index, theta_deg, phi0_mag, phi0_phase,
  phi90_mag, phi90_phase``. Angles in the CSV are degrees; magnitudes are
  linear voltage patterns relative to isotropic (|mag|^2 is element gain);
  phases are degrees.

Usage::

    python -m edgefem.contract out_dir/

writes ``golden_array_package.json`` and ``golden_patterns.csv``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

import pyedgefem

#: Bump when the fixture content changes; consumers pin this via metadata.
FIXTURE_REVISION = 1

_N_X, _N_Y = 2, 2
_SPACING_M = 0.015  # half wavelength at 10 GHz
_FREQS_HZ = (10.0e9, 10.5e9)
_THETA_DEG = np.arange(-90.0, 91.0, 5.0)


def canonical_array_package() -> pyedgefem.ArrayPackage:
    """The canonical 2x2 fixture package with hand-pinned values."""
    pkg = pyedgefem.ArrayPackage()
    pkg.model_name = f"golden-2x2-rev{FIXTURE_REVISION}"
    pkg.reference_impedance = 50.0

    positions = []
    for ix in range(_N_X):
        for iy in range(_N_Y):
            positions.append(np.array([ix * _SPACING_M, iy * _SPACING_M, 0.0]))
    pkg.element_positions = positions

    pkg.frequencies = list(_FREQS_HZ)

    # Hand-pinned S-matrices: -15 dB return on the diagonal, -20 dB adjacent
    # coupling, -25 dB diagonal coupling, with small phases so the [re, im]
    # serialization is actually exercised. Scaled slightly at the second
    # frequency so per-frequency indexing is testable.
    n = _N_X * _N_Y
    s = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            xi, yi = divmod(i, _N_Y)
            xj, yj = divmod(j, _N_Y)
            dist = math.hypot(xi - xj, yi - yj)
            if i == j:
                s[i, j] = 10 ** (-15 / 20) * np.exp(1j * math.radians(30.0))
            elif dist == 1.0:
                s[i, j] = 10 ** (-20 / 20) * np.exp(-1j * math.radians(45.0))
            else:
                s[i, j] = 10 ** (-25 / 20) * np.exp(1j * math.radians(10.0))
    pkg.S_matrices = [s, s * 0.95]

    return pkg


def canonical_patterns() -> list[pyedgefem.EmbeddedPattern]:
    """One embedded pattern per port at the first frequency.

    Cosine-shaped linear voltage magnitude with a small port-dependent phase
    ramp, so consumers can verify port indexing, the degree conversion, and
    the two-cut structure.
    """
    package = canonical_array_package()
    patterns = []
    for port, position in enumerate(package.element_positions):
        pat = pyedgefem.EmbeddedPattern()
        pat.port_index = port
        pat.frequency = _FREQS_HZ[0]
        pat.element_position = position

        phi0, phi90 = [], []
        for theta_deg in _THETA_DEG:
            theta = math.radians(float(theta_deg))
            for cuts, scale in ((phi0, 1.0), (phi90, 0.9)):
                point = pyedgefem.FFPoint2D()
                point.theta = theta  # radians, per the C++ struct
                point.magnitude = scale * max(math.cos(theta), 0.0)
                point.phase = math.radians(2.0 * port + 0.1 * float(theta_deg))
                cuts.append(point)
        pat.pattern_phi0 = phi0
        pat.pattern_phi90 = phi90
        patterns.append(pat)
    return patterns


def write_golden_array_package(out_dir: str | Path) -> tuple[Path, Path]:
    """Write the golden JSON + patterns CSV pair into *out_dir*."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / "golden_array_package.json"
    csv_path = out / "golden_patterns.csv"

    pyedgefem.export_array_package_json(str(json_path), canonical_array_package())
    pyedgefem.export_patterns_csv(str(csv_path), canonical_patterns())
    return json_path, csv_path


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "."
    written = write_golden_array_package(target)
    for path in written:
        print(f"wrote {path}")
