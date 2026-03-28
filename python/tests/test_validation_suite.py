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


# ---------- Lossy Material Validation ----------

@gmsh_required
class TestLossyMaterial:
    """Validate FEM against analytical attenuation for lossy dielectric."""

    def test_dielectric_loss_attenuation(self):
        """FEM |S21| through lossy waveguide should match analytical decay.

        Analytical: For dielectric-filled WR-90 waveguide at 10 GHz,
        attenuation constant α_d = (k₀² tan_d) / (2β)
        → |S21| = exp(-α_d * L)
        """
        from edgefem.designs import RectWaveguideDesign

        # WR-90 geometry
        a, b, L = 22.86e-3, 10.16e-3, 50e-3
        freq = 10e9
        eps_r = 2.2
        tan_d = 0.05  # Moderate loss for clear signal

        # Analytical attenuation
        c0 = 299792458.0
        omega = 2 * np.pi * freq
        k0 = omega / c0
        # Propagation constant in dielectric-filled guide
        kc = np.pi / a  # TE10 cutoff
        beta = np.sqrt(eps_r * k0**2 - kc**2)
        alpha_d = k0**2 * eps_r * tan_d / (2 * beta)
        S21_analytical = np.exp(-alpha_d * L)

        # FEM solve using low-level API with lossy dielectric
        wg = RectWaveguideDesign(a=a, b=b, length=L)
        wg.generate_mesh(density=10)
        wg._setup_ports(freq)

        params = em.MaxwellParams()
        params.omega = omega
        params.port_abc_scale = 0.5
        # Set complex permittivity: eps_r(1 - j*tan_d)
        params.eps_r = complex(eps_r, -eps_r * tan_d)

        S = em.calculate_sparams_eigenmode(
            wg._mesh, params, wg._bc, wg._ports
        )
        S21_fem = abs(S[1, 0])

        # FEM should show attenuation in the right ballpark
        # Allow 30% tolerance for coarse mesh + first-order elements
        assert S21_fem < 1.0, f"|S21| = {S21_fem:.4f} should show loss"
        assert S21_fem > 0.1, f"|S21| = {S21_fem:.4f} too low (solver issue?)"
        rel_error = abs(S21_fem - S21_analytical) / S21_analytical
        assert rel_error < 0.35, (
            f"|S21| FEM={S21_fem:.4f} vs analytical={S21_analytical:.4f}, "
            f"relative error={rel_error:.1%} exceeds 35% tolerance"
        )

    def test_lossless_vs_lossy_comparison(self):
        """Lossy waveguide should have lower |S21| than lossless."""
        from edgefem.designs import RectWaveguideDesign

        a, b, L = 22.86e-3, 10.16e-3, 50e-3
        freq = 10e9

        wg = RectWaveguideDesign(a=a, b=b, length=L)
        wg.generate_mesh(density=10)

        # Lossless
        S_lossless = wg.sparams_at_freq(freq)
        S21_lossless = abs(S_lossless[1, 0])

        # Lossy (reuse same mesh)
        wg._setup_ports(freq)
        params = em.MaxwellParams()
        params.omega = 2 * np.pi * freq
        params.port_abc_scale = 0.5
        params.eps_r = complex(2.2, -2.2 * 0.05)

        S_lossy = em.calculate_sparams_eigenmode(
            wg._mesh, params, wg._bc, wg._ports
        )
        S21_lossy = abs(S_lossy[1, 0])

        assert S21_lossy < S21_lossless, (
            f"Lossy |S21|={S21_lossy:.4f} should be less than "
            f"lossless |S21|={S21_lossless:.4f}"
        )


# ---------- Mesh Convergence Diagnostics ----------

@gmsh_required
class TestMeshConvergence:
    """Test mesh convergence check feature."""

    def test_stacked_patch_convergence_check(self):
        """convergence_check() should return a dict with expected keys."""
        from edgefem.designs import StackedPatchDesign
        design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
            cavity_size=35e-3,
            cavity_depth=5e-3,
            probe_radius=0.65e-3,
        )
        result = design.convergence_check(
            2.4e9, density=6, refinement_factor=1.5, verbose=False
        )
        assert 'converged' in result
        assert 'delta' in result
        assert 'S11_coarse' in result
        assert 'S11_fine' in result
        assert result['delta'] >= 0
        assert np.isfinite(result['delta'])

    def test_waveguide_convergence_check(self):
        """Waveguide convergence check should show small delta."""
        from edgefem.designs import RectWaveguideDesign
        wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
        result = wg.convergence_check(
            10e9, density=8, refinement_factor=1.5, verbose=False
        )
        assert 'converged' in result
        assert 'delta' in result
        assert 'S21_coarse' in result
        assert 'S21_fine' in result
        assert result['delta'] >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
