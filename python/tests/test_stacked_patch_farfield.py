"""
Far-field evaluation of StackedPatchDesign: the vectorized Stratton-Chu sum
against the loop implementation it replaced, and the per-(frequency, grid)
pattern cache.

Requires Gmsh to be installed and in PATH.
"""

import subprocess

import numpy as np
import pytest

try:
    import pyedgefem as em
    from edgefem.designs import StackedPatchDesign
    from edgefem.designs import stacked_patch
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


HAS_GMSH = _check_gmsh() if HAS_EDGEFEM else False
gmsh_required = pytest.mark.skipif(not HAS_GMSH, reason="Gmsh not available")

FREQ = 2.4e9


def reference_pattern(design, freq, n_theta, n_phi):
    """The per-direction, per-face loop that radiation_pattern() used before
    it was vectorized. Kept verbatim as the reference for the 1e-10 check.
    Returns (theta_grid, phi_grid, E_theta, E_phi)."""
    solution = design._last_solutions[freq]
    omega = 2 * np.pi * freq
    c0 = 299792458.0
    k0 = omega / c0

    huygens = em.extract_huygens_surface(
        design._mesh, solution, design._TAG_HUYGENS, omega
    )

    normals_fixed = []
    E_doubled = []
    H_zero = []
    for i in range(len(huygens.n)):
        n = huygens.n[i]
        E = huygens.E_tan[i]
        if n[2] < 0:
            normals_fixed.append(-n)
            E_doubled.append(-2.0 * E)
        else:
            normals_fixed.append(n)
            E_doubled.append(2.0 * E)
        H_zero.append(np.zeros(3, dtype=complex))

    theta_rad = np.linspace(0, np.pi / 2, n_theta).tolist()
    phi_rad = np.linspace(0, 2 * np.pi, n_phi).tolist()

    E_theta_out = np.zeros((n_theta, n_phi), dtype=complex)
    E_phi_out = np.zeros((n_theta, n_phi), dtype=complex)
    theta_grid = np.zeros((n_theta, n_phi))
    phi_grid = np.zeros((n_theta, n_phi))

    for ti, th in enumerate(theta_rad):
        sin_th = np.sin(th)
        cos_th = np.cos(th)
        obliquity = cos_th

        for pi_idx, ph in enumerate(phi_rad):
            theta_grid[ti, pi_idx] = th
            phi_grid[ti, pi_idx] = ph

            rhat = np.array([sin_th * np.cos(ph),
                             sin_th * np.sin(ph), cos_th])
            th_hat = np.array([cos_th * np.cos(ph),
                               cos_th * np.sin(ph), -sin_th])
            ph_hat = np.array([-np.sin(ph), np.cos(ph), 0.0])

            Efar = np.zeros(3, dtype=complex)
            for s in range(len(huygens.r)):
                M = -np.cross(normals_fixed[s], E_doubled[s])
                phase = np.exp(-1j * k0 * np.dot(rhat, huygens.r[s]))
                term = 1j * k0 * np.cross(np.cross(rhat, M), rhat)
                Efar += term * phase * huygens.area[s] * obliquity

            E_theta_out[ti, pi_idx] = np.dot(Efar, th_hat)
            E_phi_out[ti, pi_idx] = np.dot(Efar, ph_hat)

    return theta_grid, phi_grid, E_theta_out, E_phi_out


def _max_rel_err(a, b):
    """Largest |a - b| relative to the largest |b|."""
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


@pytest.fixture(scope="module")
def solved_design():
    """A coarse 2.4 GHz FR-4 patch, meshed and solved once for the module."""
    if not HAS_GMSH:
        pytest.skip("Gmsh not available")
    design = StackedPatchDesign(
        patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
        substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
        probe_radius=0.65e-3,
    )
    design.generate_mesh(density=6, f_min=FREQ)
    design.simulate(FREQ)
    return design


@pytest.fixture
def design(solved_design):
    """The solved design with an empty pattern cache."""
    solved_design._clear_pattern_cache()
    return solved_design


@gmsh_required
class TestVectorizedMatchesLoop:

    def test_matches_reference_loop(self, design):
        """The vectorized sum equals the loop to 1e-10 relative to the peak."""
        n_theta, n_phi = 19, 13
        pattern = design.radiation_pattern(FREQ, n_theta=n_theta, n_phi=n_phi)
        th, ph, e_th, e_ph = reference_pattern(design, FREQ, n_theta, n_phi)

        assert pattern.E_theta.shape == (n_theta, n_phi)
        np.testing.assert_array_equal(pattern.theta_grid, th)
        np.testing.assert_array_equal(pattern.phi_grid, ph)
        # Not a vacuous comparison: the pattern has a real broadside peak.
        assert np.max(np.abs(e_th)) > 0.0
        assert _max_rel_err(pattern.E_theta, e_th) < 1e-10
        assert _max_rel_err(pattern.E_phi, e_ph) < 1e-10

    def test_matches_reference_across_chunks(self, design, monkeypatch):
        """A chunk budget that splits the directions into several chunks (and
        leaves a short last chunk) gives the same field."""
        n_theta, n_phi = 11, 9
        # 99 directions; a budget for 4 directions per chunk gives 25 chunks.
        n_faces = len(em.extract_huygens_surface(
            design._mesh, design._last_solutions[FREQ],
            design._TAG_HUYGENS, 2 * np.pi * FREQ).area)
        per_dir = design._FARFIELD_BYTES_PER_ENTRY * n_faces
        monkeypatch.setattr(design, "_FARFIELD_CHUNK_BYTES", 4 * per_dir)
        pattern = design.radiation_pattern(FREQ, n_theta=n_theta, n_phi=n_phi)
        _, _, e_th, e_ph = reference_pattern(design, FREQ, n_theta, n_phi)

        assert _max_rel_err(pattern.E_theta, e_th) < 1e-10
        assert _max_rel_err(pattern.E_phi, e_ph) < 1e-10

    def test_single_direction_grid(self, design):
        """n_theta = n_phi = 1 (one direction at broadside) still works."""
        pattern = design.radiation_pattern(FREQ, n_theta=1, n_phi=1)
        _, _, e_th, e_ph = reference_pattern(design, FREQ, 1, 1)
        assert pattern.E_theta.shape == (1, 1)
        assert _max_rel_err(pattern.E_theta, e_th) < 1e-10
        assert _max_rel_err(pattern.E_phi, e_ph) < 1e-10


@gmsh_required
class TestPatternCache:

    @staticmethod
    def _count_extractions(monkeypatch):
        """Count calls into the Huygens extraction, the first step of every
        far-field evaluation."""
        calls = []
        real = stacked_patch.em.extract_huygens_surface

        def counting(*args, **kwargs):
            calls.append(args)
            return real(*args, **kwargs)

        monkeypatch.setattr(stacked_patch.em, "extract_huygens_surface",
                            counting)
        return calls

    def test_directivity_reuses_pattern(self, design, monkeypatch):
        """directivity() after radiation_pattern() on the same explicit grid
        does not evaluate the far field again."""
        calls = self._count_extractions(monkeypatch)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        assert len(calls) == 1

        D = design.directivity(FREQ, n_theta=19, n_phi=13)
        assert len(calls) == 1
        assert np.isfinite(D) and D > 1.0

    def test_directivity_default_ignores_earlier_patterns(self, design):
        """directivity() with no grid is the 91 x 73 value whether or not a
        pattern on another grid was computed first."""
        fresh = design.directivity(FREQ)
        design._clear_pattern_cache()
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        after_other_grid = design.directivity(FREQ)
        assert after_other_grid == fresh
        assert (FREQ, 91, 73) in design._pattern_cache

    def test_directivity_without_pattern_uses_default_grid(
            self, design, monkeypatch):
        calls = self._count_extractions(monkeypatch)
        design.directivity(FREQ)
        assert len(calls) == 1
        assert (FREQ, 91, 73) in design._pattern_cache

    def test_repeat_call_is_a_hit(self, design, monkeypatch):
        calls = self._count_extractions(monkeypatch)
        first = design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        second = design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        assert len(calls) == 1
        np.testing.assert_array_equal(first.E_theta, second.E_theta)
        # Callers get their own arrays: editing one result does not edit the
        # cached pattern.
        second.E_theta[:] = 0.0
        third = design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        np.testing.assert_array_equal(third.E_theta, first.E_theta)

    def test_other_grid_is_a_miss(self, design, monkeypatch):
        calls = self._count_extractions(monkeypatch)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=9)
        assert len(calls) == 2

    def test_new_solution_invalidates(self, design, monkeypatch):
        """Solving the same frequency again replaces the field, so the cached
        pattern is recomputed."""
        calls = self._count_extractions(monkeypatch)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        design.simulate(FREQ)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        assert len(calls) == 2

    def test_sweep_resolve_invalidates(self, design, monkeypatch):
        """A sweep that keeps the field at a frequency with a cached pattern
        replaces that field, so the pattern is recomputed from the new one."""
        calls = self._count_extractions(monkeypatch)
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        old_solution = design._last_solutions[FREQ]
        design.frequency_sweep(FREQ, FREQ, n_points=1, keep_fields=[FREQ])
        assert design._last_solutions[FREQ] is not old_solution
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        assert len(calls) == 2
        assert calls[1][1] is design._last_solutions[FREQ]

    def test_solutions_are_read_only(self, design):
        """_store_solution() is the only way to replace a field, so no write
        can bypass the pattern invalidation."""
        with pytest.raises(TypeError):
            design._last_solutions[FREQ] = design._last_solutions[FREQ]

    def test_new_solution_drops_its_patterns(self, design):
        """Replacing the field at a frequency frees its cached patterns."""
        design.radiation_pattern(FREQ, n_theta=19, n_phi=13)
        design.radiation_pattern(FREQ, n_theta=5, n_phi=7)
        assert len(design._pattern_cache) == 2
        design._store_solution(FREQ, design._last_solutions[FREQ])
        assert design._pattern_cache == {}

    def test_cache_holds_at_most_the_cap(self, design, monkeypatch):
        """The cache keeps the most recent _PATTERN_CACHE_MAX patterns."""
        monkeypatch.setattr(design, "_PATTERN_CACHE_MAX", 2)
        for n_phi in (3, 5, 7):
            design.radiation_pattern(FREQ, n_theta=3, n_phi=n_phi)
        assert list(design._pattern_cache) == [(FREQ, 3, 5), (FREQ, 3, 7)]

    def test_new_mesh_invalidates(self, design):
        """generate_mesh() drops cached patterns along with the solutions.
        Uses a throwaway design so the module's solved design is untouched."""
        other = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4}],
        )
        other._pattern_cache[(FREQ, 19, 13)] = object()
        other.generate_mesh(density=6)
        assert other._pattern_cache == {}
