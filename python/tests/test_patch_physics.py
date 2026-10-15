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
- test_probe_follows_probe_radius, test_domain_starts_at_ground_*,
  test_air_box_*, test_cavity_arguments_are_deprecated: feed and domain
  geometry (probe from probe_radius, no cavity, air box lambda/2).
- test_patch_known_answer_resonance (slow): resonance against the
  cavity model.
"""

import math
import warnings

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
        probe_radius=0.5e-3,
        probe_offset=(0.0, FR4_PROBE_Y_OFFSET),
    )


def _tet_volume_by_tag(mesh):
    """Total tetrahedron volume (m^3) per volume physical tag."""
    xyz = mesh.node_coords()
    conn = mesh.tet_nodes()
    phys = mesh.tet_phys()
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
# thickness. The probe post cut out of the layer removes 0.006% here
# (0.85 mm square, 1.6 mm tall), inside the 1% tolerance. On 2c4d96d the single
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
    )
    design.generate_mesh(density=6, design_freq=2.4e9,
                         output_path=str(tmp_path / "stacked.msh"))
    volumes = _tet_volume_by_tag(design.mesh)
    _assert_substrate_volumes(volumes, design.ground_plane_size,
                              [s['thickness'] for s in substrates],
                              "StackedPatchDesign (2 layers)")


# ---------------------------------------------------------------------------
# Feed and domain geometry
# ---------------------------------------------------------------------------
#
# The probe is a square post whose thin-wire equivalent radius (0.59 x
# side) equals probe_radius, at every mesh density; on 2c4d96d its size
# came from the mesh size (a 5.1 mm square at density 6, 2.5 mm at 12).
# The port strip spans the hole from the probe face to the ground on the
# +x side: (b - a) x 2a, with a, b the half-sides of probe and hole.

def _surface_area_by_tag(mesh):
    """Total triangle area (m^2) per surface physical tag."""
    xyz = mesh.node_coords()
    conn = mesh.tri_nodes()
    phys = mesh.tri_phys()
    q0, q1, q2 = (xyz[conn[:, k]] for k in range(3))
    area = np.linalg.norm(np.cross(q1 - q0, q2 - q0), axis=1) / 2.0
    return {int(tag): float(area[phys == tag].sum()) for tag in np.unique(phys)}


# The last case puts the probe 0.6 mm from the patch edge, so the hole
# (half-side 0.98 mm) crosses the edge.
NEAR_EDGE_Y_OFFSET = -FR4_PATCH['patch_length'] / 2 + 0.6e-3


@pytest.mark.parametrize("probe_radius, density, y_offset", [
    (0.5e-3, 6, FR4_PROBE_Y_OFFSET), (0.5e-3, 9, FR4_PROBE_Y_OFFSET),
    (0.8e-3, 6, FR4_PROBE_Y_OFFSET), (0.5e-3, 6, NEAR_EDGE_Y_OFFSET)])
def test_probe_follows_probe_radius(tmp_path, probe_radius, density,
                                    y_offset):
    patch = PatchAntennaDesign(**FR4_PATCH)
    patch.set_probe_feed(y_offset=y_offset, probe_radius=probe_radius)
    patch.generate_mesh(density=density, design_freq=FR4_DESIGN_FREQ,
                        output_path=str(tmp_path / "patch.msh"))
    areas = _surface_area_by_tag(patch.mesh)
    side = probe_radius / 0.59
    hole_side = side * math.exp(50.0 / 60.0)
    h = FR4_PATCH['substrate_height']
    S = StackedPatchDesign
    assert areas[S._TAG_PROBE] == pytest.approx(4 * side * h, rel=1e-6)
    assert areas[S._TAG_PORT] == pytest.approx(
        (hole_side - side) / 2 * side, rel=1e-6)
    gp = patch.ground_plane_size
    assert areas[S._TAG_GROUND] == pytest.approx(gp ** 2 - hole_side ** 2,
                                                 rel=1e-6)
    assert areas[S._TAG_PATCH_BASE] == pytest.approx(
        FR4_PATCH['patch_width'] * FR4_PATCH['patch_length'], rel=1e-6)


def test_domain_starts_at_ground_and_air_box_is_half_wavelength(tmp_path):
    # No cavity: no node below the ground plane. The air box reaches
    # lambda/2 at f_min above the patch.
    f_min = 2.0e9
    patch = _fr4_patch()
    patch.generate_mesh(density=6, design_freq=FR4_DESIGN_FREQ, f_min=f_min,
                        output_path=str(tmp_path / "patch.msh"))
    z = patch.mesh.node_coords()[:, 2]
    assert z.min() == pytest.approx(0.0, abs=1e-12)
    expected_top = FR4_PATCH['substrate_height'] + 0.5 * C0 / f_min
    assert z.max() == pytest.approx(expected_top, rel=1e-9)


def test_air_box_below_half_wavelength_warns(tmp_path):
    # Default f_min is 0.8 x design_freq: no warning at 0.85 x, one below.
    patch = _fr4_patch()
    patch.generate_mesh(density=4, design_freq=FR4_DESIGN_FREQ,
                        output_path=str(tmp_path / "patch.msh"))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        patch.simulate(0.85 * FR4_DESIGN_FREQ)
    with pytest.warns(UserWarning, match="lambda/2"):
        patch.simulate(0.75 * FR4_DESIGN_FREQ)


def test_cavity_arguments_are_deprecated():
    with pytest.warns(DeprecationWarning, match="cavity_size"):
        StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
            cavity_size=35e-3, cavity_depth=5e-3)


# ---------------------------------------------------------------------------
# Known answer: resonance of the 2.45 GHz FR-4 patch
# ---------------------------------------------------------------------------
#
# Mesh density 12, as in examples/run_patch_fullwave.py. The sweep runs
# 2.25 to 2.70 GHz in 25 MHz steps (19 points); the air box is lambda/2
# at 2.25 GHz (f_min). The cavity-model resonance
# is 2.467 GHz, so the +/-5% band (2.344 to 2.590 GHz) lies inside the
# sweep, and a minimum at either end of the sweep cannot pass.
# The result depends on the near-patch mesh size,
# StackedPatchDesign._PATCH_SIZE_X_DENSITY: 2.350 GHz (-4.7%) at 6.0 (h/2
# near the patch), 2.375 GHz (-3.7%) at 4.0 (h/3, the default since the
# convergence study in docs/validation/patch_2p45ghz.md). openEMS puts
# this patch at 2.389 GHz, 3.2% below the cavity model.
# Cost with UMFPACK: 112 s and 8.4 GB peak (process, including meshing)
# at h/3; 36 s and 6.5 GB at h/2. With the SparseLU fallback before PEC
# elimination (h/2) it took about 2 min per point, 8.3 GB and 39 min.
# Deselected by default; run with `pytest -m slow`.

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
                        design_freq=FR4_DESIGN_FREQ, f_min=f_lo,
                        output_path=str(tmp_path / "patch_known_answer.msh"))
    freqs, _, s11 = patch.frequency_sweep(f_lo, f_hi, n)

    mag_db = 20 * np.log10(np.abs(np.asarray(s11)) + 1e-30)
    i_min = int(np.argmin(mag_db))
    f_res = float(freqs[i_min])
    table = ", ".join(f"{f / 1e9:.3f}:{d:.2f}" for f, d in zip(freqs, mag_db))
    detail = (f"cavity model {f_cav / 1e9:.4f} GHz; |S11| dB by GHz: {table}")
    print(f"f_res {f_res / 1e9:.4f} GHz, min |S11| {mag_db[i_min]:.2f} dB, "
          f"{100 * (f_res - f_cav) / f_cav:+.2f}% from the {detail}")

    assert 0 < i_min < n - 1, (
        f"|S11| minimum at the sweep edge ({f_res / 1e9:.3f} GHz); {detail}")
    assert mag_db[i_min] < MIN_S11_DB, (
        f"min |S11| = {mag_db[i_min]:.2f} dB, need < {MIN_S11_DB} dB; {detail}")
    assert abs(f_res - f_cav) / f_cav < F_RES_RTOL, (
        f"f_res = {f_res / 1e9:.4f} GHz is "
        f"{100 * (f_res - f_cav) / f_cav:+.1f}% from the cavity model; {detail}")
