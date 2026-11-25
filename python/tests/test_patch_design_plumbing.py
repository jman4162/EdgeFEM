"""
PatchAntennaDesign plumbing: edge feeds in FEM, load_mesh, the z0 argument
on the FEM path, the analytical impedance model and pattern export.
"""

import csv
import subprocess
import warnings

import numpy as np
import pytest

try:
    from edgefem.designs import PatchAntennaDesign, StackedPatchDesign
    from edgefem.designs.stacked_patch import _AperturePattern
    HAS_EDGEFEM = True
except ImportError:
    HAS_EDGEFEM = False

pytestmark = pytest.mark.skipif(not HAS_EDGEFEM, reason="pyedgefem not built")


def _check_gmsh():
    try:
        subprocess.run(["gmsh", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


gmsh_required = pytest.mark.skipif(
    not (HAS_EDGEFEM and _check_gmsh()), reason="Gmsh not available")

FREQ = 2.4e9


def _fr4_patch():
    return PatchAntennaDesign(patch_length=29.5e-3, patch_width=38e-3,
                              substrate_height=1.6e-3, substrate_eps_r=4.4)


def test_edge_feed_is_refused_by_fem():
    """The FEM backend models probe feeds only. An edge feed used to become
    a probe at the patch centre (the TM010 null) without a word."""
    patch = _fr4_patch()
    patch.set_edge_feed(inset_length=5e-3, line_width=3e-3)
    with pytest.raises(NotImplementedError, match="edge feed"):
        patch.generate_mesh(density=6)


class TestAnalyticalImpedance:

    def test_impedance_is_continuous(self):
        """No jump anywhere in 0.6 to 1.4 f_r (the old model set X = 0
        outside +-20 %)."""
        patch = _fr4_patch()
        patch.set_probe_feed()
        fr = patch._estimated_fr
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            Z = np.array([patch.input_impedance(f)
                          for f in np.linspace(0.6 * fr, 1.4 * fr, 2001)])
        assert np.max(np.abs(np.diff(Z))) < 0.05 * np.max(np.abs(Z))

    def test_resonance_and_bandwidth(self):
        """Real Z at f_r, capacitive above (parallel RLC), and a VSWR-2
        bandwidth of a few percent for 1.6 mm FR-4 (the old model gave
        about 40 %)."""
        patch = _fr4_patch()
        patch.set_probe_feed()
        fr = patch._estimated_fr
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            Z_r = patch.input_impedance(fr)
            Z_hi = patch.input_impedance(1.02 * fr)
            _, _, bw = patch.bandwidth()
        assert abs(Z_r.imag) < 1e-9 * abs(Z_r)
        assert Z_hi.imag < 0.0
        assert 0.5 < bw < 10.0


def test_analytical_pattern_exports_csv(tmp_path):
    patch = _fr4_patch()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        path = tmp_path / "pattern.csv"
        patch.export_pattern_csv(str(path), 2.45e9)
    with open(path) as f:
        rows = list(csv.reader(f))
    assert rows[0][:2] == ["theta_deg", "phi_deg"]
    assert len(rows) == 1 + 91 * 73


def test_aperture_pattern_has_pattern_db():
    th, ph = np.meshgrid(np.linspace(0, np.pi / 2, 5),
                         np.linspace(0, 2 * np.pi, 7), indexing="ij")
    pat = _AperturePattern(th, ph, np.cos(th).astype(complex),
                           np.zeros_like(th, dtype=complex))
    dB = pat.pattern_dB()
    assert dB.shape == th.shape
    assert dB.max() == 0.0


@gmsh_required
class TestFemPlumbing:

    @pytest.fixture(scope="class")
    def solved(self, tmp_path_factory):
        path = str(tmp_path_factory.mktemp("mesh") / "patch.msh")
        patch = _fr4_patch()
        patch.set_probe_feed()
        patch.generate_mesh(density=6, output_path=path, f_min=FREQ)
        return patch, path

    def test_load_mesh_enables_fem(self, solved):
        _, path = solved
        other = _fr4_patch()
        other.set_probe_feed()
        other.load_mesh(path)
        assert other.has_fem
        S11, _ = other.simulate(FREQ)
        assert np.isfinite(S11)

    def test_return_loss_honours_z0(self, solved):
        patch, _ = solved
        Z = patch.input_impedance(FREQ)
        for z0 in (50.0, 75.0):
            want = 20 * np.log10(abs((Z - z0) / (Z + z0)))
            assert abs(patch.return_loss(FREQ, z0=z0) - want) < 1e-9

    def test_sweep_honours_z0(self, solved):
        patch, _ = solved
        f, Z, S11 = patch.frequency_sweep(FREQ, FREQ, n_points=1, z0=75.0)
        assert abs(S11[0] - (Z[0] - 75.0) / (Z[0] + 75.0)) < 1e-9
