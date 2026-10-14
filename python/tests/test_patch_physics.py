"""Physics checks for the probe-fed patch designs.

These tests come from section D of the EdgeFEM investigation of
2026-10-12. On commit 2c4d96d every patch run returned S11 = -0.99x at
every frequency because three defects stack: the design discards the
ports that normalize_port_weights returns, the C++ lumped-port term acts
as a reactance, and the geometry has no substrate volume and shorts the
probe through a PEC cavity. The earlier patch tests only checked
|S11| <= 1.05, which that output satisfies.

- test_stacked_patch_s11_is_not_a_short_*: observable effect of the
  non-physical port term (S11 pinned at -1). Since the gap-port fix the
  designs no longer call normalize_port_weights.
- test_*_substrate_volume*: the substrate tag must fill the substrate.
- test_patch_known_answer_resonance (slow): resonance against the
  cavity model.
"""

import math

import numpy as np
import pytest

try:
    import pyedgefem as em  # noqa: F401
    from edgefem.designs import PatchAntennaDesign, StackedPatchDesign
    HAS_EDGEFEM = True
except ImportError:
    HAS_EDGEFEM = False

try:
    import gmsh  # noqa: F401
    HAS_GMSH = True
except ImportError:
    HAS_GMSH = False

pytestmark = [
    pytest.mark.skipif(not HAS_EDGEFEM, reason="pyedgefem not built"),
    pytest.mark.skipif(not HAS_GMSH, reason="gmsh Python module not installed"),
]

C0 = 299792458.0

# 2.45 GHz probe-fed patch on FR-4, from examples/run_patch_fullwave.py.
FR4_PATCH = dict(
    patch_length=28.6e-3,
    patch_width=37.3e-3,
    substrate_height=1.6e-3,
    substrate_eps_r=4.4,
    substrate_tan_d=0.02,
)
FR4_PROBE_Y_OFFSET = -9.0e-3
FR4_DESIGN_FREQ = 2.45e9


def _fr4_patch():
    patch = PatchAntennaDesign(**FR4_PATCH)
    patch.set_probe_feed(y_offset=FR4_PROBE_Y_OFFSET)
    return patch


def _fr4_stacked_design():
    """The FR-4 patch built directly as a StackedPatchDesign.

    Same arguments PatchAntennaDesign passes to its backend.
    """
    p = FR4_PATCH
    return StackedPatchDesign(
        patches=[{
            'width': p['patch_width'],
            'length': p['patch_length'],
            'height': p['substrate_height'],
        }],
        substrates=[{
            'thickness': p['substrate_height'],
            'eps_r': p['substrate_eps_r'],
            'tan_d': p['substrate_tan_d'],
        }],
        cavity_size=max(p['patch_width'], p['patch_length']) + 10e-3,
        cavity_depth=5e-3,
        probe_radius=0.5e-3,
        probe_offset=(0.0, FR4_PROBE_Y_OFFSET),
    )


def _tet_volume_by_tag(mesh):
    """Total tetrahedron volume (m^3) per volume physical tag."""
    # Each attribute access copies the whole vector; read each once.
    nodes = mesh.nodes
    node_index = mesh.nodeIndex
    tets = mesh.tets
    xyz = np.array([n.xyz for n in nodes], dtype=float)
    conn = np.array([[node_index[c] for c in t.conn[:4]] for t in tets])
    phys = np.array([t.phys for t in tets])
    p0, p1, p2, p3 = (xyz[conn[:, k]] for k in range(4))
    vol = np.abs(np.einsum(
        'ij,ij->i', p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0
    return {int(tag): float(vol[phys == tag].sum()) for tag in np.unique(phys)}


def _cavity_model_resonance(length, width, h, eps_r):
    """TM010 resonance of a rectangular patch (cavity model).

    Balanis, Antenna Theory (3rd ed.), Sec. 14.2, with Hammerstad's
    fringing length extension:
        eps_reff = (eps_r + 1)/2 + (eps_r - 1)/2 * (1 + 12 h / W)^(-1/2)
        dL = 0.412 h (eps_reff + 0.3)(W/h + 0.264)
             / ((eps_reff - 0.258)(W/h + 0.8))
        f_r = c0 / (2 (L + 2 dL) sqrt(eps_reff))
    For the FR-4 patch (L 28.6 mm, W 37.3 mm, h 1.6 mm, eps_r 4.4):
    eps_reff = 4.08, dL = 0.74 mm, f_r = 2.467 GHz.
    """
    eps_reff = (eps_r + 1) / 2 + (eps_r - 1) / 2 / math.sqrt(1 + 12 * h / width)
    dl = (0.412 * h * (eps_reff + 0.3) * (width / h + 0.264)
          / ((eps_reff - 0.258) * (width / h + 0.8)))
    return C0 / (2 * (length + 2 * dl) * math.sqrt(eps_reff))


# ---------------------------------------------------------------------------
# Normalize round trip (observable behavior)
# ---------------------------------------------------------------------------
#
# With the normalized weights dropped and the reactive port term, the
# shipped path gives S11 = -0.9947, -0.9938, -0.9921 at 2.2, 2.45 and
# 2.7 GHz (density 6). A design whose port reaches the solver correctly
# cannot sit at S11 = -1 at every frequency: that is a short circuit at
# the port. The investigation's physical readout on the same (still
# flawed) geometry gave S11 = -0.71+0.69j, -0.56+0.81j, -0.29+0.93j.

S11_SHORT_TOL = 0.01
ROUND_TRIP_FREQS = [2.2e9, 2.45e9, 2.7e9]


def _assert_not_pinned_at_short(s11, label):
    s11 = np.asarray(s11, dtype=complex)
    dist = np.abs(s11 + 1.0)
    assert not np.all(dist < S11_SHORT_TOL), (
        f"{label}: S11 within {S11_SHORT_TOL} of -1 at every frequency "
        f"(S11 = {np.round(s11, 4).tolist()}); the port does not reach "
        f"the solver as a physical load"
    )


def test_stacked_patch_s11_is_not_a_short_simulate(tmp_path):
    design = _fr4_stacked_design()
    design.generate_mesh(density=6, design_freq=FR4_DESIGN_FREQ,
                         output_path=str(tmp_path / "patch.msh"))
    s11 = [design.simulate(f)[0] for f in ROUND_TRIP_FREQS]
    _assert_not_pinned_at_short(s11, "StackedPatchDesign.simulate")


def test_stacked_patch_s11_is_not_a_short_sweep(tmp_path):
    design = _fr4_stacked_design()
    design.generate_mesh(density=6, design_freq=FR4_DESIGN_FREQ,
                         output_path=str(tmp_path / "patch.msh"))
    _, _, s11 = design.frequency_sweep(
        ROUND_TRIP_FREQS[0], ROUND_TRIP_FREQS[-1], len(ROUND_TRIP_FREQS))
    _assert_not_pinned_at_short(s11, "StackedPatchDesign.frequency_sweep")


# ---------------------------------------------------------------------------
# Substrate coverage
# ---------------------------------------------------------------------------
#
# The substrate under the ground plane footprint must be meshed and tagged
# as substrate: total tet volume per substrate tag = ground area x layer
# thickness. A probe cut out of the layer removes at most 0.3% here (the
# 5 mm square probe that 2c4d96d builds at density 6; 0.01% for a 0.5 mm
# radius probe), inside the 1% tolerance. On 2c4d96d the single
# patch's substrate tag holds only the 67 tets inside the probe
# (investigation A3) and the second layer of the dual stack has no tets.

SUBSTRATE_VOLUME_RTOL = 0.01


def _assert_substrate_volumes(volumes, ground_size, thicknesses, label):
    base = StackedPatchDesign._TAG_SUBSTRATE_BASE
    for i, t in enumerate(thicknesses):
        expected = ground_size ** 2 * t
        got = volumes.get(base + i, 0.0)
        assert got == pytest.approx(expected, rel=SUBSTRATE_VOLUME_RTOL), (
            f"{label}: substrate layer {i} (tag {base + i}) volume "
            f"{got * 1e9:.1f} mm^3, expected ground area x thickness = "
            f"{expected * 1e9:.1f} mm^3; volumes by tag (mm^3): "
            f"{ {k: round(v * 1e9, 1) for k, v in volumes.items()} }"
        )


def test_patch_antenna_substrate_volume(tmp_path):
    patch = _fr4_patch()
    patch.generate_mesh(density=6, design_freq=FR4_DESIGN_FREQ,
                        output_path=str(tmp_path / "patch.msh"))
    volumes = _tet_volume_by_tag(patch.mesh)
    _assert_substrate_volumes(volumes, patch.ground_plane_size,
                              [FR4_PATCH['substrate_height']],
                              "PatchAntennaDesign")


def test_stacked_patch_substrate_volumes(tmp_path):
    substrates = [
        {'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02},
        {'thickness': 3.2e-3, 'eps_r': 1.05},
    ]
    design = StackedPatchDesign(
        patches=[
            {'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3},
            {'width': 32e-3, 'length': 27e-3, 'height': 4.8e-3},
        ],
        substrates=substrates,
        cavity_size=35e-3,
        cavity_depth=5e-3,
    )
    design.generate_mesh(density=6, design_freq=2.4e9,
                         output_path=str(tmp_path / "stacked.msh"))
    volumes = _tet_volume_by_tag(design.mesh)
    _assert_substrate_volumes(volumes, design.ground_plane_size,
                              [s['thickness'] for s in substrates],
                              "StackedPatchDesign (2 layers)")


# ---------------------------------------------------------------------------
# Known answer: resonance of the 2.45 GHz FR-4 patch
# ---------------------------------------------------------------------------
#
# Mesh density 12, as in examples/run_patch_fullwave.py. The sweep runs
# 2.25 to 2.70 GHz in 25 MHz steps (19 points). The cavity-model resonance
# is 2.467 GHz, so the +/-5% band (2.344 to 2.590 GHz) lies inside the
# sweep, and a minimum at either end of the sweep cannot pass.
# Cost: one LU per point at density 12 (about 35 s and 6 GB each on an
# idle machine), so this runs for minutes. Deselected by default; run with
# `pytest -m slow`.

KNOWN_ANSWER_DENSITY = 12
KNOWN_ANSWER_SWEEP = (2.25e9, 2.70e9, 19)
F_RES_RTOL = 0.05
MIN_S11_DB = -6.0


@pytest.mark.slow
def test_patch_known_answer_resonance(tmp_path):
    f_cav = _cavity_model_resonance(
        FR4_PATCH['patch_length'], FR4_PATCH['patch_width'],
        FR4_PATCH['substrate_height'], FR4_PATCH['substrate_eps_r'])
    f_lo, f_hi, n = KNOWN_ANSWER_SWEEP
    assert f_lo < (1 - F_RES_RTOL) * f_cav and f_hi > (1 + F_RES_RTOL) * f_cav

    patch = _fr4_patch()
    patch.generate_mesh(density=KNOWN_ANSWER_DENSITY,
                        design_freq=FR4_DESIGN_FREQ,
                        output_path=str(tmp_path / "patch_known_answer.msh"))
    freqs, _, s11 = patch.frequency_sweep(f_lo, f_hi, n)

    mag_db = 20 * np.log10(np.abs(np.asarray(s11)) + 1e-30)
    i_min = int(np.argmin(mag_db))
    f_res = float(freqs[i_min])
    table = ", ".join(f"{f / 1e9:.3f}:{d:.2f}" for f, d in zip(freqs, mag_db))
    detail = (f"cavity model {f_cav / 1e9:.4f} GHz; |S11| dB by GHz: {table}")

    assert 0 < i_min < n - 1, (
        f"|S11| minimum at the sweep edge ({f_res / 1e9:.3f} GHz); {detail}")
    assert mag_db[i_min] < MIN_S11_DB, (
        f"min |S11| = {mag_db[i_min]:.2f} dB, need < {MIN_S11_DB} dB; {detail}")
    assert abs(f_res - f_cav) / f_cav < F_RES_RTOL, (
        f"f_res = {f_res / 1e9:.4f} GHz is "
        f"{100 * (f_res - f_cav) / f_cav:+.1f}% from the cavity model; {detail}")
