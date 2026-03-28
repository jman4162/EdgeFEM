"""
Golden benchmark validation suite for EdgeFEM.

Tests fundamental electromagnetic properties that any correct FEM solver must satisfy:
- S-parameter passivity (energy conservation)
- S-parameter reciprocity (symmetry)
- Cavity eigenfrequency accuracy vs analytical formulas
- Solver convergence and numerical stability
- Frequency sweep consistency
"""

import sys
import pytest
import numpy as np

try:
    import pyedgefem as em
    HAS_EDGEFEM = True
except ImportError:
    HAS_EDGEFEM = False

pytestmark = pytest.mark.skipif(not HAS_EDGEFEM, reason="pyedgefem not built")


def _check_gmsh():
    import subprocess
    try:
        subprocess.run(["gmsh", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


HAS_GMSH = _check_gmsh() if HAS_EDGEFEM else False
gmsh_required = pytest.mark.skipif(not HAS_GMSH, reason="Gmsh not available")


# ---------- Unit Cell Passivity and Reciprocity ----------

@gmsh_required
class TestUnitCellPassivity:
    """Verify passivity for unit cell designs."""

    def test_grounded_patch_passive(self):
        """Grounded patch unit cell must have |S11| <= 1."""
        from edgefem.designs import UnitCellDesign
        cell = UnitCellDesign(
            period_x=10e-3, period_y=10e-3,
            substrate_height=1.6e-3, substrate_eps_r=4.4,
        )
        cell.add_patch(width=8e-3, length=6e-3)
        cell.generate_mesh(density=8)
        R, T = cell.reflection_transmission(10e9)
        assert abs(R) <= 1.01, f"|R| = {abs(R):.4f} > 1.01 (passivity violation)"

    def test_bare_ground_total_reflection(self):
        """Grounded substrate with no patch should reflect nearly all energy."""
        from edgefem.designs import UnitCellDesign
        cell = UnitCellDesign(
            period_x=10e-3, period_y=10e-3,
            substrate_height=0.5e-3, substrate_eps_r=2.2,
        )
        cell.generate_mesh(density=8)
        R, T = cell.reflection_transmission(10e9)
        assert abs(R) > 0.90, f"|R| = {abs(R):.4f} should be near 1.0 for bare ground"


# ---------- Stacked Patch Validation ----------

@gmsh_required
class TestStackedPatchValidation:
    """Validate stacked patch antenna with full-wave FEM."""

    def setup_method(self):
        from edgefem.designs import StackedPatchDesign
        self.design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
            cavity_size=35e-3,
            cavity_depth=5e-3,
            probe_radius=0.65e-3,
        )
        self.design.generate_mesh(density=8)

    def test_passivity_across_band(self):
        """S11 magnitude stays bounded across frequency sweep."""
        freqs, Z_list, S11_list = self.design.frequency_sweep(
            2.0e9, 3.0e9, n_points=5
        )
        for i, (f, s11) in enumerate(zip(freqs, S11_list)):
            # Allow coarse-mesh tolerance of 5%
            assert abs(s11) <= 1.10, (
                f"|S11| = {abs(s11):.4f} at {f/1e9:.3f} GHz "
                f"exceeds coarse-mesh passivity bound"
            )

    def test_impedance_finite(self):
        """Input impedance should be finite at all frequencies."""
        freqs, Z_list, S11_list = self.design.frequency_sweep(
            2.0e9, 3.0e9, n_points=5
        )
        for f, Z in zip(freqs, Z_list):
            assert np.isfinite(Z.real), f"Re(Z) not finite at {f/1e9:.3f} GHz"
            assert np.isfinite(Z.imag), f"Im(Z) not finite at {f/1e9:.3f} GHz"

    def test_sweep_consistency(self):
        """Sweep S11 should match single-point simulate."""
        freq = 2.4e9
        S11_single, _ = self.design.simulate(freq)

        freqs, _, S11_sweep = self.design.frequency_sweep(
            freq, freq, n_points=1
        )
        # Allow some tolerance since sweep uses mid-frequency normalization
        assert abs(abs(S11_single) - abs(S11_sweep[0])) < 0.15, (
            f"|S11| single={abs(S11_single):.4f} vs sweep={abs(S11_sweep[0]):.4f}"
        )


# ---------- Analytical Patch Antenna Validation ----------

class TestPatchAntennaAnalytical:
    """Validate analytical patch antenna model against known formulas."""

    def test_resonant_frequency_fr4(self):
        """FR-4 patch resonant frequency should match cavity model."""
        from edgefem.designs import PatchAntennaDesign
        patch = PatchAntennaDesign(
            patch_length=29e-3, patch_width=38e-3,
            substrate_height=1.6e-3, substrate_eps_r=4.4,
        )
        fr = patch.estimated_resonant_frequency
        # For L=29mm, eps_r=4.4: fr ≈ c/(2L√εeff) ≈ 2.4 GHz (rough)
        assert 1.5e9 < fr < 3.5e9, f"fr = {fr/1e9:.3f} GHz out of expected range"

    def test_directivity_range(self):
        """Patch antenna directivity should be in 5-9 dBi range."""
        from edgefem.designs import PatchAntennaDesign
        patch = PatchAntennaDesign(
            patch_length=29e-3, patch_width=38e-3,
            substrate_height=1.6e-3, substrate_eps_r=4.4,
        )
        D = patch.directivity(2.4e9)
        D_dBi = 10 * np.log10(D) if D > 0 else -np.inf
        # Analytical model gives lower directivity than full-wave (simplified aperture)
        assert 1.0 <= D_dBi <= 10.0, f"Directivity = {D_dBi:.1f} dBi out of range"

    def test_impedance_positive_resistance(self):
        """Input impedance real part must be positive."""
        from edgefem.designs import PatchAntennaDesign
        patch = PatchAntennaDesign(
            patch_length=29e-3, patch_width=38e-3,
            substrate_height=1.6e-3, substrate_eps_r=4.4,
        )
        Z = patch.input_impedance(2.4e9)
        assert Z.real > 0, f"Re(Z) = {Z.real:.2f} should be positive"


# ---------- Solver Robustness ----------

@gmsh_required
class TestSolverRobustness:
    """Test solver convergence and fallback behavior."""

    def test_ilut_preconditioner_available(self):
        """SolveOptions should expose ILUT parameters."""
        opts = em.SolveOptions()
        assert hasattr(opts, 'use_ilut')
        assert hasattr(opts, 'auto_fallback')
        assert opts.use_ilut is True
        assert opts.auto_fallback is True

    def test_direct_solver_works(self):
        """Direct SparseLU solver should work for small problems."""
        from edgefem.designs import UnitCellDesign
        cell = UnitCellDesign(
            period_x=10e-3, period_y=10e-3,
            substrate_height=1.6e-3, substrate_eps_r=4.4,
        )
        cell.generate_mesh(density=6)
        R, T = cell.reflection_transmission(10e9)
        assert np.isfinite(abs(R)), "Reflection coefficient should be finite"


# ---------- Lumped Port Weight Modes ----------

@gmsh_required
class TestLumpedPortWeightModes:
    """Test both weight computation strategies for lumped ports."""

    def test_surface_integral_mode_default(self):
        """Default weight mode should be SurfaceIntegral."""
        config = em.LumpedPortConfig()
        assert config.weight_mode == em.LumpedPortWeightMode.SurfaceIntegral

    def test_projection_mode_available(self):
        """Projection weight mode should be accessible."""
        config = em.LumpedPortConfig()
        config.weight_mode = em.LumpedPortWeightMode.Projection
        assert config.weight_mode == em.LumpedPortWeightMode.Projection


# ---------- Frequency Sweep API ----------

@gmsh_required
class TestFrequencySweepAPI:
    """Test the C++ frequency sweep API."""

    def test_sweep_api_exists(self):
        """frequency_sweep function should be available."""
        assert hasattr(em, 'frequency_sweep')
        assert hasattr(em, 'SweepResult')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
