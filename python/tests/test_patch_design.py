"""
Integration tests for PatchAntennaDesign.

Tests the design class workflow and validates that computed antenna parameters
are physically reasonable. Uses the analytical model path (no FEM solve
required for impedance), so these tests run fast.
"""

import numpy as np
import pytest

pyedgefem = pytest.importorskip("pyedgefem", reason="pyedgefem not built")
from edgefem.designs.patch_antenna import PatchAntennaDesign


@pytest.fixture
def patch_2ghz():
    """2.4 GHz patch antenna on FR-4 substrate."""
    return PatchAntennaDesign(
        patch_length=29e-3,
        patch_width=38e-3,
        substrate_height=1.6e-3,
        substrate_eps_r=4.4,
    )


class TestPatchBasic:
    def test_estimated_resonance(self, patch_2ghz):
        """Estimated resonant frequency should be near 2.4 GHz for these dimensions."""
        fr = patch_2ghz.estimated_resonant_frequency
        # Transmission line model gives approximate value; allow 20% tolerance
        assert 1.8e9 < fr < 3.2e9, f"Estimated fr = {fr/1e9:.2f} GHz, expected ~2.4 GHz"

    def test_repr(self, patch_2ghz):
        r = repr(patch_2ghz)
        assert "PatchAntennaDesign" in r
        assert "mm" in r


class TestPatchImpedance:
    def test_impedance_real_positive(self, patch_2ghz):
        """Input impedance real part must be positive (passive device)."""
        patch_2ghz.set_probe_feed()
        Z_in = patch_2ghz.input_impedance(2.4e9)
        assert Z_in.real > 0, f"Re(Z_in) = {Z_in.real:.1f}, must be positive"

    def test_impedance_reasonable_magnitude(self, patch_2ghz):
        """Input impedance should be in a reasonable range (1-1000 ohms)."""
        patch_2ghz.set_probe_feed()
        Z_in = patch_2ghz.input_impedance(2.4e9)
        assert 1 < abs(Z_in) < 1000, f"|Z_in| = {abs(Z_in):.1f}, out of range"

    def test_impedance_near_resonance_low_reactance(self, patch_2ghz):
        """At resonance, reactance should be small relative to resistance."""
        patch_2ghz.set_probe_feed()
        fr = patch_2ghz.estimated_resonant_frequency
        Z_in = patch_2ghz.input_impedance(fr)
        # At resonance, |X| should be small compared to R
        assert abs(Z_in.imag) < abs(Z_in.real) * 2, (
            f"Z_in = {Z_in:.1f} — reactance too large at resonance"
        )


class TestPatchReturnLoss:
    def test_return_loss_is_negative_db(self, patch_2ghz):
        """Return loss should be negative (in dB)."""
        patch_2ghz.set_probe_feed()
        rl = patch_2ghz.return_loss(2.4e9)
        assert rl < 0, f"Return loss = {rl:.1f} dB, should be negative"

    def test_return_loss_bounded(self, patch_2ghz):
        """Return loss should be between -60 dB and 0 dB."""
        patch_2ghz.set_probe_feed()
        rl = patch_2ghz.return_loss(2.4e9)
        assert -60 < rl < 0, f"Return loss = {rl:.1f} dB, out of range"


class TestPatchRadiation:
    def test_directivity_positive(self, patch_2ghz):
        """Directivity should be > 1 (> 0 dBi) for a patch antenna."""
        D = patch_2ghz.directivity(2.4e9)
        D_dBi = 10 * np.log10(D)
        assert D_dBi > 0, f"Directivity = {D_dBi:.1f} dBi, expected > 0"

    def test_directivity_reasonable(self, patch_2ghz):
        """Patch directivity is typically 5-9 dBi."""
        D = patch_2ghz.directivity(2.4e9)
        D_dBi = 10 * np.log10(D)
        assert 2 < D_dBi < 15, f"Directivity = {D_dBi:.1f} dBi, expected 5-9 dBi"

    def test_pattern_shape(self, patch_2ghz):
        """Pattern should have maximum at broadside (theta=0)."""
        pattern = patch_2ghz.radiation_pattern(2.4e9, n_theta=19, n_phi=9)
        pwr = pattern.power_pattern()
        # Broadside is at theta=0 (first row)
        max_idx = np.unravel_index(np.argmax(pwr), pwr.shape)
        assert max_idx[0] == 0, "Pattern maximum not at broadside"


class TestPatchFeedTypes:
    def test_probe_feed(self, patch_2ghz):
        patch_2ghz.set_probe_feed(x_offset=0, y_offset=-5e-3)
        Z_in = patch_2ghz.input_impedance(2.4e9)
        assert Z_in.real > 0

    def test_edge_feed(self, patch_2ghz):
        patch_2ghz.set_edge_feed(inset_length=5e-3)
        Z_in = patch_2ghz.input_impedance(2.4e9)
        assert Z_in.real > 0


def _check_gmsh():
    import subprocess
    try:
        subprocess.run(["gmsh", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


HAS_GMSH = _check_gmsh()
gmsh_required = pytest.mark.skipif(not HAS_GMSH, reason="Gmsh not available")


@gmsh_required
class TestPatchFEM:
    """Tests using full-wave FEM simulation."""

    def setup_method(self):
        self.patch = PatchAntennaDesign(
            patch_length=29e-3,
            patch_width=38e-3,
            substrate_height=1.6e-3,
            substrate_eps_r=4.4,
            substrate_tan_d=0.02,
        )

    def test_fem_mesh_generation(self):
        """FEM mesh generation via StackedPatchDesign backend works."""
        self.patch.generate_mesh(density=8)
        assert self.patch.mesh is not None
        assert self.patch.has_fem

    def test_fem_simulate(self):
        """FEM simulate returns S11 with reasonable magnitude."""
        self.patch.generate_mesh(density=8)
        S11, sol = self.patch.simulate(2.4e9)
        assert np.isfinite(abs(S11))
        # Coarse mesh tolerance
        assert abs(S11) <= 1.10, f"|S11| = {abs(S11):.4f} exceeds tolerance"

    def test_fem_input_impedance(self):
        """FEM input impedance is finite."""
        self.patch.generate_mesh(density=8)
        Z_in = self.patch.input_impedance(2.4e9)
        assert np.isfinite(Z_in.real)
        assert np.isfinite(Z_in.imag)

    def test_fem_return_loss(self):
        """FEM return loss is finite."""
        self.patch.generate_mesh(density=8)
        rl = self.patch.return_loss(2.4e9)
        assert np.isfinite(rl)

    def test_fem_radiation_pattern(self):
        """FEM radiation pattern has broadside maximum."""
        self.patch.generate_mesh(density=8)
        pattern = self.patch.radiation_pattern(2.4e9, n_theta=19, n_phi=9)
        power = pattern.power_pattern()
        assert power.shape[0] > 0
        # Maximum should be near broadside (theta=0)
        max_idx = np.unravel_index(np.argmax(power), power.shape)
        theta_max = pattern.theta_grid[max_idx]
        assert theta_max < np.pi / 4, (
            f"Pattern max at theta={np.degrees(theta_max):.1f} deg"
        )

    def test_analytical_fallback_without_mesh(self):
        """Without mesh, input_impedance uses analytical path."""
        # No mesh generated — should use analytical
        Z_in = self.patch.input_impedance(2.4e9)
        assert Z_in.real > 0  # Analytical should work without mesh
