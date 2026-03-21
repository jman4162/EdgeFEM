"""
Tests for StackedPatchDesign.

Requires Gmsh to be installed and in PATH.
"""

import sys
import pytest
import numpy as np

try:
    import pyedgefem as em
    from edgefem.designs import StackedPatchDesign
    HAS_EDGEFEM = True
except ImportError:
    HAS_EDGEFEM = False

pytestmark = pytest.mark.skipif(not HAS_EDGEFEM, reason="pyedgefem not built")


def _check_gmsh():
    """Check if Gmsh is available."""
    import subprocess
    try:
        subprocess.run(["gmsh", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


HAS_GMSH = _check_gmsh() if HAS_EDGEFEM else False
gmsh_required = pytest.mark.skipif(not HAS_GMSH, reason="Gmsh not available")


@gmsh_required
class TestSinglePatch:
    """Tests for single-patch cavity-fed antenna."""

    def setup_method(self):
        self.design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
            cavity_size=35e-3,
            cavity_depth=5e-3,
            probe_radius=0.65e-3,
        )

    def test_mesh_generation(self):
        """Mesh generation produces valid mesh."""
        self.design.generate_mesh(density=8)
        mesh = self.design.mesh
        assert mesh is not None
        assert mesh.num_tets() > 0
        assert mesh.num_edges() > 0

    def test_simulate_s11(self):
        """S11 magnitude is less than 1 (passive device)."""
        self.design.generate_mesh(density=8)
        S11, solution = self.design.simulate(2.4e9)
        assert abs(S11) <= 1.05, f"|S11| = {abs(S11):.4f} > 1.05 (passivity violation)"
        assert len(solution) > 0

    def test_input_impedance(self):
        """Input impedance has positive real part."""
        self.design.generate_mesh(density=8)
        Z_in = self.design.input_impedance(2.4e9)
        assert Z_in.real > 0, f"Re(Z_in) = {Z_in.real:.2f} should be positive"

    def test_return_loss(self):
        """Return loss is finite and negative."""
        self.design.generate_mesh(density=8)
        rl = self.design.return_loss(2.4e9)
        assert np.isfinite(rl)
        assert rl <= 0, f"Return loss = {rl:.2f} dB should be <= 0"


@gmsh_required
class TestDualPatch:
    """Tests for dual-patch stacked antenna."""

    def setup_method(self):
        self.design = StackedPatchDesign(
            patches=[
                {'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3},
                {'width': 32e-3, 'length': 27e-3, 'height': 4.8e-3},
            ],
            substrates=[
                {'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02},
                {'thickness': 3.2e-3, 'eps_r': 1.05},
            ],
            cavity_size=35e-3,
            cavity_depth=5e-3,
        )

    def test_mesh_generation(self):
        """Dual-patch mesh generation works."""
        self.design.generate_mesh(density=8)
        assert self.design.mesh is not None
        assert self.design.mesh.num_tets() > 0

    def test_simulate(self):
        """Dual-patch simulation completes."""
        self.design.generate_mesh(density=8)
        S11, solution = self.design.simulate(2.4e9)
        assert abs(S11) <= 1.05
        assert np.isfinite(abs(S11))


@gmsh_required
class TestRadiationPattern:
    """Tests for radiation pattern computation."""

    @pytest.mark.xfail(reason="Huygens surface (tag 60) not yet generated in mesh")
    def test_broadside_maximum(self):
        """Radiation pattern has maximum near broadside (theta=0)."""
        design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
            cavity_size=35e-3,
            cavity_depth=5e-3,
        )
        design.generate_mesh(density=8)
        pattern = design.radiation_pattern(2.4e9, n_theta=37, n_phi=19)

        # Pattern should have maximum near theta=0 (broadside)
        power = pattern.power_pattern()
        assert power.shape[0] > 0
        assert power.shape[1] > 0

        # Find theta index of maximum
        max_idx = np.unravel_index(np.argmax(power), power.shape)
        theta_max = pattern.theta_grid[max_idx]
        # Broadside is theta=0; allow up to 45 degrees due to coarse mesh
        assert theta_max < np.pi / 4, (
            f"Pattern maximum at theta={np.degrees(theta_max):.1f} deg, "
            f"expected near broadside"
        )

    @pytest.mark.xfail(reason="Huygens surface (tag 60) not yet generated in mesh")
    def test_directivity_range(self):
        """Directivity is in expected range for patch antenna (5-10 dBi)."""
        design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
            cavity_size=35e-3,
            cavity_depth=5e-3,
        )
        design.generate_mesh(density=8)
        D = design.directivity(2.4e9)
        D_dBi = 10 * np.log10(D)
        # Patch antennas typically have 5-10 dBi directivity
        # Allow wider range for coarse mesh: 2-15 dBi
        assert 2 <= D_dBi <= 15, f"Directivity = {D_dBi:.1f} dBi out of range"


@gmsh_required
def test_repr():
    """Test string representation."""
    design = StackedPatchDesign(
        patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
        substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
        cavity_size=35e-3,
        cavity_depth=5e-3,
    )
    s = repr(design)
    assert "StackedPatchDesign" in s
    assert "1 patch" in s


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
