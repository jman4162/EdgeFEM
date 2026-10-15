"""
Integration tests for UnitCellDesign.

Tests periodic unit cell simulation workflow. Uses coarse meshes for speed.
Validates that reflection/transmission coefficients are physical.
"""

from pathlib import Path

import numpy as np
import pytest

pyedgefem = pytest.importorskip("pyedgefem", reason="pyedgefem not built")
from edgefem.designs.unit_cell import UnitCellDesign


@pytest.fixture
def simple_cell():
    """Simple unit cell with patch on grounded substrate."""
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1.0e-3,
        substrate_eps_r=3.5,
    )
    cell.add_patch(width=4e-3, length=4e-3)
    cell.generate_mesh(density=8)
    return cell


@pytest.fixture
def bare_cell():
    """Bare unit cell (no patch) — should have high reflection from ground plane."""
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1.0e-3,
        substrate_eps_r=3.5,
    )
    cell.generate_mesh(density=8)
    return cell


class TestUnitCellBasic:
    def test_mesh_generated(self, simple_cell):
        assert simple_cell.mesh is not None

    def test_not_fss(self, simple_cell):
        """Default config has ground plane — not an FSS."""
        assert not simple_cell.is_fss

    def test_repr(self, simple_cell):
        r = repr(simple_cell)
        assert "UnitCellDesign" in r
        assert "period" in r


class TestUnitCellReflection:
    def test_reflection_magnitude_bounded(self, simple_cell):
        """Reflection coefficient must satisfy |R| <= 1 (passive)."""
        R, T = simple_cell.reflection_transmission(10e9)
        assert abs(R) <= 1.05, f"|R| = {abs(R):.4f}, violates passivity"

    def test_grounded_no_transmission(self, simple_cell):
        """Grounded cell should have zero transmission."""
        R, T = simple_cell.reflection_transmission(10e9)
        assert abs(T) < 0.05, f"|T| = {abs(T):.4f}, expected ~0 for grounded cell"

    def test_bare_ground_high_reflection(self, bare_cell):
        """Bare grounded substrate should reflect most energy."""
        R, T = bare_cell.reflection_transmission(10e9)
        assert abs(R) > 0.5, f"|R| = {abs(R):.4f}, expected high for bare ground"

    def test_reflection_returns_complex(self, simple_cell):
        R, T = simple_cell.reflection_transmission(10e9)
        assert isinstance(R, complex)
        assert isinstance(T, (complex, float, int))


class TestUnitCellSurfaceImpedance:
    def test_surface_impedance_finite(self, simple_cell):
        """Surface impedance should be finite for a loaded cell."""
        Zs = simple_cell.surface_impedance(10e9)
        assert np.isfinite(Zs.real), f"Re(Zs) = {Zs.real}, not finite"
        assert np.isfinite(Zs.imag), f"Im(Zs) = {Zs.imag}, not finite"

    def test_surface_impedance_positive_real(self, simple_cell):
        """Surface impedance real part should be positive (lossy or lossless)."""
        Zs = simple_cell.surface_impedance(10e9)
        # Real part should be non-negative for passive structure
        assert Zs.real >= -1.0, f"Re(Zs) = {Zs.real:.1f}, expected >= 0"


class TestUnitCellFSS:
    def test_fss_creation(self):
        """FSS has air below and no ground plane."""
        cell = UnitCellDesign(
            period_x=5e-3,
            period_y=5e-3,
            substrate_height=0.5e-3,
            substrate_eps_r=3.5,
            air_height_below=15e-3,
            has_ground_plane=False,
        )
        assert cell.is_fss

    def test_effective_parameters_requires_fss(self, simple_cell):
        """effective_parameters should raise for non-FSS cell."""
        with pytest.raises(ValueError, match="FSS"):
            simple_cell.effective_parameters(10e9)


# ---------------------------------------------------------------------------
# Floquet ports keep the 1.0.0 port terms
# ---------------------------------------------------------------------------
# build_wave_port gives Modal ports. calculate_sparams_periodic uses the
# 1.0.0 modal terms for them (load w w^H / Z0, source 2 w / sqrt(Z0)); only
# build_lumped_port ports get the gap-port terms. These S-parameters are not
# validated, and every call emits a RuntimeWarning.
#
# data/unit_cell_bare_d8.msh: UnitCellDesign(period_x=5e-3, period_y=5e-3,
# substrate_height=1e-3, substrate_eps_r=3.5) without a patch,
# generate_mesh(density=8), gmsh 4.15.0, 466 tets. The air above the slab is
# the default lambda/2 at 10 GHz, 14.990 mm, so the port plane is at
# z = 15.990 mm.
_BARE_MESH = str(Path(__file__).parent / "data" / "unit_cell_bare_d8.msh")

# R = S11 from EdgeFEM 1.0.0 (Python 3.12 wheel) on that mesh, 2026-10-15.
_R_1P0 = {
    8e9: complex(-0.9906263162642917, -0.13659978596117378),
    10e9: complex(-0.9963762268305602, -0.08505536201111075),
    12e9: complex(-0.9983423078818191, -0.057555506193612625),
}


def _bare_cell_from_mesh():
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1.0e-3,
        substrate_eps_r=3.5,
    )
    cell.load_mesh(_BARE_MESH)
    return cell


def grounded_slab_reflection(freq, eps_r, h, d):
    """Normal-incidence reflection of a PEC-backed dielectric slab.

    Time convention exp(+j omega t). The slab (thickness h, relative
    permittivity eps_r) on a PEC ground has input impedance
    Z_s = j eta_d tan(k_d h), eta_d = eta0 / sqrt(eps_r),
    k_d = k0 sqrt(eps_r). At the slab face
    Gamma_s = (Z_s - eta0) / (Z_s + eta0) = -(1 - j t) / (1 + j t),
    t = tan(k_d h) / sqrt(eps_r), so |Gamma_s| = 1 and
    arg Gamma_s = pi - 2 atan(t). Referred to a port plane d above the slab
    face, Gamma = Gamma_s exp(-2 j k0 d).
    """
    c0 = 299792458.0
    k0 = 2 * np.pi * freq / c0
    t = np.tan(k0 * np.sqrt(eps_r) * h) / np.sqrt(eps_r)
    gamma_s = -(1 - 1j * t) / (1 + 1j * t)
    return complex(gamma_s * np.exp(-2j * k0 * d))


class TestFloquetPort1p0:
    def test_runtime_warning(self):
        cell = _bare_cell_from_mesh()
        with pytest.warns(RuntimeWarning, match="not validated"):
            cell.reflection_transmission(10e9)

    @pytest.mark.filterwarnings("ignore:S-parameters computed with modal")
    def test_matches_1p0(self):
        """Branch R equals the 1.0.0 R on the same mesh to 1e-9."""
        cell = _bare_cell_from_mesh()
        for f, r_ref in _R_1P0.items():
            R, _ = cell.reflection_transmission(f)
            assert abs(R - r_ref) < 1e-9, f"{f/1e9} GHz: R = {R}, 1.0.0 {r_ref}"


@pytest.mark.known_issue
@pytest.mark.xfail(
    strict=True,
    reason="Floquet port S-parameters are not validated (1.0.0 port terms)",
)
@pytest.mark.filterwarnings("ignore:S-parameters computed with modal")
def test_grounded_slab_reflection_phase():
    """Known answer: bare grounded slab under a Floquet port, normal incidence.

    eps_r = 3.5, h = 1 mm, port plane d = 14.990 mm above the slab
    (grounded_slab_reflection). Analytic arg(R): -127.67 deg at 8 GHz,
    +155.07 deg at 10 GHz (2 k0 d = 360 deg there), +77.58 deg at 12 GHz;
    |R| = 1.

    Measured on 2026-10-15 with this mesh (density 8): -172.15, -175.12 and
    -176.70 deg, |R| = 1.0000 (the 1.0.0 values; the gap-port terms that
    briefly applied to Floquet ports gave R = +0.9995 at 0 deg at all three
    frequencies). The phase errors are 44.5, 29.8 and 105.7 deg, so the
    test fails; it is marked xfail(strict=True), and a fix of the Floquet
    port turns it into an XPASS failure that asks for the marks to be
    removed. The 25 deg tolerance allows for this coarse mesh; R = +1 at
    0 deg misses by 127.7, 155.1 and 77.6 deg.
    """
    cell = _bare_cell_from_mesh()
    d = cell.air_height_above
    errors = []
    for f in (8e9, 10e9, 12e9):
        R, _ = cell.reflection_transmission(f)
        ref = grounded_slab_reflection(f, 3.5, 1.0e-3, d)
        err = np.degrees(np.angle(R / ref))
        errors.append(abs(err))
        assert abs(abs(R) - 1.0) < 0.05, f"{f/1e9} GHz: |R| = {abs(R):.4f}"
    assert max(errors) < 25.0, f"phase errors (deg): {errors}"
