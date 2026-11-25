"""
Far-field conventions of the Python aperture model and the analytical patch
pattern, checked on synthetic apertures with closed-form answers. The
solver's time convention is exp(+j omega t).

For an aperture field E_a on z = 0 over an infinite PEC ground (n = +z),
image theory gives the magnetic current M = -2 n x E_a radiating in free
space, and (Balanis 12-10, r-normalized, exp(-j k r)/r removed)
    E_theta = -(j k / 4 pi) L_phi,  E_phi = (j k / 4 pi) L_theta,
    L = sum M area exp(+j k rhat . r').
"""

import numpy as np
import pytest

try:
    from edgefem.designs import StackedPatchDesign
    from edgefem.designs.stacked_patch import _AperturePattern
    from edgefem.designs import PatchAntennaDesign
    HAS_EDGEFEM = True
except ImportError:
    HAS_EDGEFEM = False

pytestmark = pytest.mark.skipif(not HAS_EDGEFEM, reason="pyedgefem not built")

K0 = 2 * np.pi  # wavelength 1 m


def _aperture(side, n_cells, field):
    """Square aperture on z = 0 centred at the origin. field(x, y) returns
    the tangential E (3-vector) at a cell centre."""
    h = side / n_cells
    c = -side / 2 + (np.arange(n_cells) + 0.5) * h
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.stack([X.ravel(), Y.ravel(), np.zeros(X.size)], axis=1)
    E = np.array([field(x, y) for x, y in r[:, :2]], dtype=complex)
    n = np.tile([0.0, 0.0, 1.0], (r.shape[0], 1))
    area = np.full(r.shape[0], h * h)
    return r, n, E, area


def test_uniform_aperture_broadside_field_and_polarization():
    """A uniform x-polarized aperture radiates an x-polarized broadside
    field of magnitude k E0 A / (2 pi) (E0 = 1)."""
    side = 2.0
    r, n, E, area = _aperture(side, 20, lambda x, y: [1.0, 0.0, 0.0])
    E_th, E_ph = StackedPatchDesign._aperture_far_field(
        r, n, E, area, K0, np.array([0.0]), np.array([0.0, np.pi / 2]))
    expected = K0 * side * side / (2 * np.pi)
    # At theta = 0: phi = 0 has theta_hat = x, phi = 90 deg has phi_hat = -x.
    assert abs(abs(E_th[0, 0]) - expected) < 1e-9 * expected
    assert abs(E_ph[0, 0]) < 1e-9 * expected
    assert abs(abs(E_ph[0, 1]) - expected) < 1e-9 * expected
    assert abs(E_th[0, 1]) < 1e-9 * expected
    # Phase: E_theta(0, 0) = -(j k / 4 pi) L_phi with L = -2 y A.
    assert abs(E_th[0, 0] - 1j * expected) < 1e-9 * expected


def test_tilted_aperture_peaks_at_plus_theta0():
    """E_a = x exp(-j k sin(th0) x) is a wave leaving toward +x (exp(+j
    omega t)); the beam peaks at (th0, phi = 0), not phi = pi."""
    th0 = np.deg2rad(30.0)
    r, n, E, area = _aperture(
        6.0, 60, lambda x, y: [np.exp(-1j * K0 * np.sin(th0) * x), 0, 0])
    theta = np.deg2rad(np.arange(0, 91))
    E_th, E_ph = StackedPatchDesign._aperture_far_field(
        r, n, E, area, K0, theta, np.array([0.0, np.pi]))
    pwr = np.abs(E_th) ** 2 + np.abs(E_ph) ** 2
    i, j = np.unravel_index(np.argmax(pwr), pwr.shape)
    assert j == 0
    assert abs(theta[i] - th0) <= np.deg2rad(2.0)
    assert pwr[30, 1] < 1e-2 * pwr[i, j]


def test_directivity_of_cos_squared_hemisphere_is_six():
    """U = cos^2(theta) on the upper hemisphere, zero below:
    P = 2 pi / 3, D = 4 pi / P = 6."""
    theta = np.linspace(0, np.pi / 2, 181)
    phi = np.linspace(0, 2 * np.pi, 145)
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    pat = _AperturePattern(TH, PH, np.cos(TH).astype(complex),
                           np.zeros_like(TH, dtype=complex))
    D = StackedPatchDesign._compute_directivity_py(pat)
    assert abs(D - 6.0) < 0.01 * 6.0


def test_large_uniform_aperture_directivity():
    """A uniform 4 x 4 wavelength aperture: D close to 4 pi A / lambda^2."""
    side = 4.0
    r, n, E, area = _aperture(side, 40, lambda x, y: [1.0, 0.0, 0.0])
    theta = np.linspace(0, np.pi / 2, 361)
    phi = np.linspace(0, 2 * np.pi, 145)
    E_th, E_ph = StackedPatchDesign._aperture_far_field(
        r, n, E, area, K0, theta, phi)
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    D = StackedPatchDesign._compute_directivity_py(
        _AperturePattern(TH, PH, E_th, E_ph))
    D_ideal = 4 * np.pi * side * side
    assert abs(D - D_ideal) < 0.05 * D_ideal


def test_aperture_normals_from_extractor_must_point_up():
    """The model takes n = +z; a face with n = -z would flip its current."""
    r, n, E, area = _aperture(1.0, 4, lambda x, y: [1.0, 0.0, 0.0])
    n[0] = [0.0, 0.0, -1.0]
    with pytest.raises(ValueError):
        StackedPatchDesign._aperture_far_field(
            r, n, E, area, K0, np.array([0.0]), np.array([0.0]))


class TestAnalyticalPatchPattern:

    @staticmethod
    def _patch():
        p = PatchAntennaDesign(patch_length=29.5e-3, patch_width=38e-3,
                               substrate_height=1.6e-3, substrate_eps_r=4.4)
        p.set_probe_feed()
        return p

    def test_no_radiation_below_ground(self):
        with pytest.warns(UserWarning):
            pat = self._patch().radiation_pattern(2.45e9, n_theta=37, n_phi=25)
        below = pat.theta_grid > np.pi / 2 + 1e-12
        assert below.any()
        assert np.all(pat.power_pattern()[below] == 0.0)

    def test_broadside_polarized_along_resonant_length(self):
        """The TM010 field across the radiating slots points along the patch
        length (y): at broadside, E_phi(phi = 0) carries it, E_theta = 0."""
        with pytest.warns(UserWarning):
            pat = self._patch().radiation_pattern(2.45e9, n_theta=37, n_phi=25)
        assert abs(pat.E_theta[0, 0]) < 1e-9 * abs(pat.E_phi[0, 0])
        assert abs(pat.E_phi[0, 0]) > 0.0

    def test_directivity_in_textbook_range(self):
        """Two-slot model over an infinite ground: 6.0 dBi for this FR-4
        patch (5.5 to 8 dBi accepted; Balanis 14.2 gives 6 to 7 dBi for
        patches near lambda/2). With the back lobe it read 3.0 dBi."""
        with pytest.warns(UserWarning):
            D = self._patch().directivity(2.45e9)
        assert 5.5 <= 10 * np.log10(D) <= 8.0
