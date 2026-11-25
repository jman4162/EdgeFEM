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
def bare_cell():
    """Bare grounded substrate (no patch), experimental results enabled."""
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1.0e-3,
        substrate_eps_r=3.5,
        experimental=True,
    )
    cell.generate_mesh(density=8)
    return cell


class TestUnitCellBasic:
    def test_mesh_generated(self, bare_cell):
        assert bare_cell.mesh is not None

    def test_not_fss(self, bare_cell):
        """Default config has ground plane — not an FSS."""
        assert not bare_cell.is_fss

    def test_repr(self, bare_cell):
        r = repr(bare_cell)
        assert "UnitCellDesign" in r
        assert "period" in r


class TestQuarantine:
    """Without a validated Floquet port the results are refused unless the
    caller opts in, and geometry that generate_mesh() ignores is refused."""

    def test_results_need_experimental(self):
        cell = UnitCellDesign(period_x=5e-3, period_y=5e-3)
        with pytest.raises(NotImplementedError, match="experimental=True"):
            cell.reflection_transmission(10e9)
        with pytest.raises(NotImplementedError, match="experimental=True"):
            cell.surface_impedance(10e9)

    @pytest.mark.parametrize("method,args", [
        ("add_patch", (4e-3, 4e-3)),
        ("add_slot", (1e-3, 4e-3)),
        ("add_via", (0.2e-3,)),
    ])
    def test_unmeshed_geometry_is_refused(self, method, args):
        cell = UnitCellDesign(period_x=5e-3, period_y=5e-3, experimental=True)
        with pytest.raises(NotImplementedError, match="not meshed"):
            getattr(cell, method)(*args)

    def test_oblique_y_incidence_is_refused(self, bare_cell):
        """Only the x pairs are applied, so k_y != 0 cannot be honoured."""
        with pytest.raises(ValueError, match="phi"):
            bare_cell.reflection_transmission(10e9, theta=30.0, phi=90.0)


class TestUnitCellReflection:
    @pytest.mark.filterwarnings("ignore:S-parameters computed with modal")
    def test_bare_ground_reflection_magnitude(self, bare_cell):
        """Lossless grounded slab: |R| = 1 (the phase is wrong, see the
        known-answer test below)."""
        R, T = bare_cell.reflection_transmission(10e9)
        assert isinstance(R, complex)
        assert abs(abs(R) - 1.0) < 0.01, f"|R| = {abs(R):.4f}"


class TestNRW:
    """effective_parameters' extraction on analytic slab S-parameters
    (exp(+j omega t)): Gamma = (z - 1) / (z + 1), P = exp(-j k0 n d),
    S11 = Gamma (1 - P^2) / (1 - Gamma^2 P^2),
    S21 = P (1 - Gamma^2) / (1 - Gamma^2 P^2)."""

    @staticmethod
    def _slab(eps, mu, k0, d):
        n = np.sqrt(complex(eps * mu))
        if n.imag > 0:  # passive under exp(+j omega t): Im(n) <= 0
            n = -n
        z = np.sqrt(mu / eps)
        g = (z - 1) / (z + 1)
        P = np.exp(-1j * k0 * n * d)
        den = 1 - g**2 * P**2
        return g * (1 - P**2) / den, P * (1 - g**2) / den

    @pytest.mark.parametrize("eps,mu", [
        (4.0, 1.0),
        (3.5 - 0.07j, 1.0),
        (2.2, 1.5 - 0.05j),
    ])
    def test_recovers_slab(self, eps, mu):
        k0 = 2 * np.pi * 10e9 / 299792458.0
        d = 1e-3
        S11, S21 = self._slab(eps, mu, k0, d)
        e, m = UnitCellDesign._nrw(S11, S21, k0, d)
        assert abs(e - eps) < 1e-9 * abs(eps)
        assert abs(m - mu) < 1e-9 * abs(mu)

    def test_deembeds_reference_planes(self):
        """S at port planes d_top above and d_bottom below the slab faces:
        S11_port = S11 exp(-2 j k0 d_top), S21_port = S21 exp(-j k0
        (d_top + d_bottom)); _nrw_at_ports undoes that."""
        k0 = 2 * np.pi * 10e9 / 299792458.0
        d, dt, db = 1e-3, 7e-3, 5e-3
        S11, S21 = self._slab(4.0, 1.0, k0, d)
        S11p = S11 * np.exp(-2j * k0 * dt)
        S21p = S21 * np.exp(-1j * k0 * (dt + db))
        e, m = UnitCellDesign._nrw(S11p, S21p, k0, d, d_top=dt, d_bottom=db)
        assert abs(e - 4.0) < 1e-9 * 4.0
        assert abs(m - 1.0) < 1e-9


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

    def test_effective_parameters_requires_fss(self, bare_cell):
        """effective_parameters should raise for non-FSS cell."""
        with pytest.raises(ValueError, match="FSS"):
            bare_cell.effective_parameters(10e9)


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

def _bare_cell_from_mesh():
    cell = UnitCellDesign(
        period_x=5e-3,
        period_y=5e-3,
        substrate_height=1.0e-3,
        substrate_eps_r=3.5,
        experimental=True,
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

    Measured with this mesh (density 8), after the periodic-pair
    orientation fix (2026-11): -170.59, -174.06 and -175.97 deg,
    |R| = 1.0000; phase errors 42.9, 30.9 and 106.4 deg. (1.0.0, with
    wrongly oriented pairs: -172.15, -175.12, -176.70 deg.) The port terms,
    not the periodic BC, cause the error. The test is xfail(strict=True); a
    fix of the Floquet port turns it into an XPASS failure that asks for
    the marks to be removed. The 25 deg tolerance allows for this coarse
    mesh; R = +1 at 0 deg misses by 127.7, 155.1 and 77.6 deg.
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
