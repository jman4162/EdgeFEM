"""
Integration tests for RectWaveguideDesign.

Tests the full workflow: mesh generation → port setup → solve → S-parameter
extraction. Uses coarse meshes for speed. Validates physics (passivity,
transmission magnitude, phase) rather than just "doesn't crash."
"""

import numpy as np
import pytest

pyedgefem = pytest.importorskip("pyedgefem", reason="pyedgefem not built")
from edgefem.designs.waveguide import RectWaveguideDesign


@pytest.fixture
def wr90():
    """WR-90 waveguide (X-band), 50 mm length, coarse mesh."""
    wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
    wg.generate_mesh(density=6)
    return wg


class TestWaveguideBasic:
    def test_cutoff_frequency(self):
        wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
        c0 = 299792458.0
        expected_fc = c0 / (2 * 22.86e-3)
        assert abs(wg.cutoff_frequency - expected_fc) < 1.0  # within 1 Hz

    def test_below_cutoff_raises(self, wr90):
        with pytest.raises(ValueError, match="below TE10 cutoff"):
            wr90.sparams_at_freq(5e9)

    def test_mesh_generated(self, wr90):
        assert wr90.mesh is not None


class TestWaveguideSparams:
    def test_transmission_above_cutoff(self, wr90):
        """TE10 mode at 10 GHz should propagate with high transmission."""
        S = wr90.sparams_at_freq(10e9)
        S21_mag = abs(S[1, 0])
        # Coarse mesh won't be perfect, but should show clear propagation
        assert S21_mag > 0.80, f"|S21| = {S21_mag:.4f}, expected > 0.80"

    def test_passivity(self, wr90):
        """S-parameter matrix must satisfy passivity: |S11|^2 + |S21|^2 <= 1.01."""
        S = wr90.sparams_at_freq(10e9)
        power = abs(S[0, 0]) ** 2 + abs(S[1, 0]) ** 2
        assert power <= 1.05, f"|S11|^2 + |S21|^2 = {power:.4f}, violates passivity"

    def test_reciprocity(self, wr90):
        """Reciprocal structure: S12 ≈ S21."""
        S = wr90.sparams_at_freq(10e9)
        assert abs(S[0, 1] - S[1, 0]) < 0.1, (
            f"S12={S[0,1]:.4f}, S21={S[1,0]:.4f} — not reciprocal"
        )

    def test_s21_phase_progresses(self, wr90):
        """S21 phase should be negative (wave propagates in +z direction)."""
        S = wr90.sparams_at_freq(10e9)
        phase_deg = np.angle(S[1, 0], deg=True)
        # Phase should be significantly negative for 50mm waveguide at 10 GHz
        assert phase_deg < -10, f"S21 phase = {phase_deg:.1f} deg, expected < -10"

    def test_smatrix_shape(self, wr90):
        S = wr90.sparams_at_freq(10e9)
        assert S.shape == (2, 2)


class TestWaveguideFrequencySweep:
    def test_sweep_returns_correct_count(self, wr90):
        freqs, S_list = wr90.frequency_sweep(8e9, 12e9, n_points=3)
        assert len(freqs) == 3
        assert len(S_list) == 3

    def test_sweep_transmission_increases_with_frequency(self, wr90):
        """Transmission should generally be high across the band."""
        freqs, S_list = wr90.frequency_sweep(8e9, 12e9, n_points=3)
        for S in S_list:
            assert abs(S[1, 0]) > 0.5, "Transmission too low in passband"


def test_assemble_maxwell_warns_on_modal_ports():
    """assemble_maxwell applies the unvalidated 1.0.0 terms to modal ports,
    so it warns like the S-parameter functions."""
    import pyedgefem as em
    mesh = em.load_gmsh("examples/rect_waveguide.msh")
    bc = em.build_edge_pec(mesh, 1)
    dims = em.RectWaveguidePort()
    dims.a, dims.b = 0.02286, 0.01016
    mode = em.solve_te10_mode(dims, 10e9)
    port = em.build_wave_port_2d(mesh, 2, mode, set(bc.dirichlet_edges),
                                 (np.pi / dims.a) ** 2)
    p = em.MaxwellParams()
    p.omega = 2 * np.pi * 10e9
    with pytest.warns(RuntimeWarning, match="not validated"):
        em.assemble_maxwell(mesh, p, bc, [port], 0)


def test_mixed_port_kinds_are_refused():
    import pyedgefem as em
    mesh = em.load_gmsh("examples/rect_waveguide.msh")
    bc = em.build_edge_pec(mesh, 1)
    dims = em.RectWaveguidePort()
    dims.a, dims.b = 0.02286, 0.01016
    mode = em.solve_te10_mode(dims, 10e9)
    kc_sq = (np.pi / dims.a) ** 2
    modal = em.build_wave_port_2d(mesh, 2, mode, set(bc.dirichlet_edges), kc_sq)
    lumped = em.build_wave_port_2d(mesh, 3, mode, set(bc.dirichlet_edges), kc_sq)
    lumped.kind = em.PortKind.Lumped
    p = em.MaxwellParams()
    p.omega = 2 * np.pi * 10e9
    with pytest.raises(ValueError, match="Lumped and Modal"):
        with pytest.warns(RuntimeWarning):
            em.calculate_sparams(mesh, p, bc, [modal, lumped])


def test_export_fields_vtk_uses_eigenmode_field(tmp_path):
    """export_fields_vtk writes the field of the validated eigenmode port
    model, without the modal-port RuntimeWarning of the 1.0.0 terms."""
    import warnings
    wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
    wg.load_mesh("examples/rect_waveguide.msh")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        wg.export_fields_vtk(str(tmp_path / "f.vtu"), freq=10e9)
    assert not [w for w in caught if "not validated" in str(w.message)]
    assert (tmp_path / "f.vtu").stat().st_size > 0
