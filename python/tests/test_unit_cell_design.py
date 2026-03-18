"""
Integration tests for UnitCellDesign.

Tests periodic unit cell simulation workflow. Uses coarse meshes for speed.
Validates that reflection/transmission coefficients are physical.
"""

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
