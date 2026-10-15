"""
Solve counts of StackedPatchDesign: simulate() factors once per frequency,
and frequency_sweep(keep_fields=...) keeps field solutions so a later
radiation_pattern() or directivity() at a kept frequency needs no solve.
Also the numpy mesh accessors and near_field().

Requires Gmsh to be installed and in PATH.
"""

import pickle
import subprocess

import numpy as np
import pytest

try:
    import pyedgefem as em
    from edgefem.designs import PatchAntennaDesign, StackedPatchDesign
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

F_START, F_STOP, N_POINTS = 2.3e9, 2.5e9, 3


def _lu_count():
    return em._direct_factorization_count()


@pytest.fixture
def design():
    """A coarse FR-4 patch, freshly meshed (no stored solutions)."""
    if not HAS_GMSH:
        pytest.skip("Gmsh not available")
    d = StackedPatchDesign(
        patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
        substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
        probe_radius=0.65e-3,
    )
    d.generate_mesh(density=6, design_freq=2.4e9, f_min=F_START)
    return d


def _max_rel_err(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


def _old_simulate_field(design, freq):
    """The field the former simulate() computed after calculate_sparams():
    assemble_maxwell with port 0 active, then solve_linear."""
    params = em.MaxwellParams()
    params.omega = 2 * np.pi * freq
    params.use_abc = True
    params.abc_surface_tags = {design._TAG_ABC}
    sub = design.substrates[0]
    params.set_eps_r_region(design._TAG_SUBSTRATE_BASE,
                            complex(sub['eps_r'], -sub['eps_r'] * sub['tan_d']))
    opts = em.SolveOptions()
    opts.use_direct = True
    asm = em.assemble_maxwell(design._mesh, params, design._bc,
                              design._ports, 0)
    return em.solve_linear(asm.A, asm.b, opts).x.to_numpy()


@gmsh_required
class TestSimulate:

    def test_simulate_factors_once_and_keeps_the_field(self, design):
        """One factorization gives S11 and the field; the field is the one
        the former separate solve gave (same matrix, port 0 excited)."""
        before = _lu_count()
        S11, x = design.simulate(2.4e9)
        assert _lu_count() - before == 1
        assert np.isfinite(S11) and abs(S11) < 1.0
        assert design._last_solutions[2.4e9] is x
        ref = _old_simulate_field(design, 2.4e9)
        assert _max_rel_err(x.to_numpy(), ref) < 1e-10


@gmsh_required
class TestSweepKeepsFields:

    def test_default_keeps_no_fields(self, design):
        before = _lu_count()
        result = design.frequency_sweep(F_START, F_STOP, n_points=2)
        assert _lu_count() - before == 2
        assert len(design._last_solutions) == 0
        assert result.min_s11_index is None
        assert result.min_s11_frequency is None

    def test_min_s11_pattern_and_directivity_need_no_solve(self, design):
        before = _lu_count()
        result = design.frequency_sweep(
            F_START, F_STOP, n_points=N_POINTS, keep_fields="min_s11")
        freqs, _, s11 = result
        assert _lu_count() - before == N_POINTS
        f_res = float(freqs[int(np.argmin(np.abs(s11)))])
        assert list(design._last_solutions) == [f_res]
        # The sweep reports which frequency it kept.
        assert result.min_s11_index == int(np.argmin(np.abs(s11)))
        assert result.min_s11_frequency == f_res
        copy = pickle.loads(pickle.dumps(result))
        assert copy.min_s11_frequency == f_res
        np.testing.assert_array_equal(copy[2], s11)

        before = _lu_count()
        pattern = design.radiation_pattern(f_res, n_theta=10, n_phi=13)
        D = design.directivity(f_res, n_theta=10, n_phi=13)
        assert _lu_count() - before == 0
        assert np.isfinite(D) and D > 1.0

        # The kept field is the f_res, port-0 field: simulate() at f_res
        # (assembled directly instead of from K and M) gives the same
        # pattern to round-off.
        design.simulate(f_res)
        fresh = design.radiation_pattern(f_res, n_theta=10, n_phi=13)
        assert _max_rel_err(pattern.E_theta, fresh.E_theta) < 1e-8
        assert _max_rel_err(pattern.E_phi, fresh.E_phi) < 1e-8

    def test_keep_list_of_frequencies_and_min_s11(self, design):
        before = _lu_count()
        freqs, _, s11 = design.frequency_sweep(
            F_START, F_STOP, n_points=N_POINTS,
            keep_fields=[F_START, F_STOP, "min_s11"])
        assert _lu_count() - before == N_POINTS
        np.testing.assert_array_equal(
            freqs, np.linspace(F_START, F_STOP, N_POINTS))
        f_res = float(freqs[int(np.argmin(np.abs(s11)))])
        assert set(design._last_solutions) == {F_START, F_STOP, f_res}

        before = _lu_count()
        design.radiation_pattern(F_STOP, n_theta=4, n_phi=5)
        assert _lu_count() - before == 0

    def test_frequency_off_the_grid_is_refused_before_solving(self, design):
        before = _lu_count()
        with pytest.raises(ValueError, match="not a sweep frequency"):
            design.frequency_sweep(F_START, F_STOP, n_points=N_POINTS,
                                   keep_fields=[2.45e9])
        with pytest.raises(ValueError, match="not a sweep frequency"):
            design.frequency_sweep(F_START, F_STOP, n_points=N_POINTS,
                                   keep_fields=2.45e9)
        with pytest.raises(ValueError, match="keep_fields"):
            design.frequency_sweep(F_START, F_STOP, n_points=N_POINTS,
                                   keep_fields="argmin")
        assert _lu_count() - before == 0

    def test_patch_antenna_design_passes_keep_fields(self):
        if not HAS_GMSH:
            pytest.skip("Gmsh not available")
        patch = PatchAntennaDesign(
            patch_length=25e-3, patch_width=30e-3, substrate_height=1.6e-3,
            substrate_eps_r=4.4, substrate_tan_d=0.02)
        patch.generate_mesh(density=6, design_freq=2.4e9, f_min=F_START)
        freqs, _, _ = patch.frequency_sweep(
            F_START, F_START, n_points=1, keep_fields="min_s11")
        before = _lu_count()
        patch.radiation_pattern(float(freqs[0]), n_theta=4, n_phi=5)
        assert _lu_count() - before == 0


@gmsh_required
class TestDirectBackends:

    def test_umfpack_and_sparselu_agree(self, design):
        """Both direct backends give the same S11 and field on the patch
        (Dirichlet edges eliminated in both), or UmfPack is refused in a
        build without SuiteSparse."""
        design._setup_simulation(2.4e9)
        params = em.MaxwellParams()
        params.omega = 2 * np.pi * 2.4e9
        params.use_abc = True
        params.abc_surface_tags = {design._TAG_ABC}
        sub = design.substrates[0]
        params.set_eps_r_region(
            design._TAG_SUBSTRATE_BASE,
            complex(sub['eps_r'], -sub['eps_r'] * sub['tan_d']))

        def solve(backend):
            opts = em.SolveOptions()
            opts.use_direct = True
            opts.direct_backend = backend
            S, x = em.calculate_sparams_with_fields(
                design._mesh, params, design._bc, design._ports, opts)
            return S, x[0].to_numpy()

        S_lu, x_lu = solve(em.DirectBackend.SparseLU)
        if not em.umfpack_available():
            with pytest.raises(ValueError, match="SuiteSparse"):
                solve(em.DirectBackend.UmfPack)
            return
        S_umf, x_umf = solve(em.DirectBackend.UmfPack)
        assert abs(S_umf[0, 0] - S_lu[0, 0]) / abs(S_lu[0, 0]) < 1e-9
        assert _max_rel_err(x_umf, x_lu) < 1e-9
        pec = np.array(sorted(design._bc.dirichlet_edges))
        assert np.all(x_umf[pec] == 0.0) and np.all(x_lu[pec] == 0.0)


@gmsh_required
class TestMeshArrays:

    def test_arrays_match_the_element_lists(self, design):
        mesh = design.mesh
        nodes, tets, tris = mesh.nodes, mesh.tets, mesh.tris
        node_index = mesh.nodeIndex

        xyz = mesh.node_coords()
        assert xyz.shape == (len(nodes), 3) and xyz.dtype == np.float64
        np.testing.assert_array_equal(xyz, [n.xyz for n in nodes])

        np.testing.assert_array_equal(
            mesh.tet_nodes(),
            [[node_index[c] for c in t.conn[:4]] for t in tets])
        np.testing.assert_array_equal(mesh.tet_edges(),
                                      [list(t.edges) for t in tets])
        np.testing.assert_array_equal(mesh.tet_edge_orient(),
                                      [list(t.edge_orient) for t in tets])
        np.testing.assert_array_equal(mesh.tet_phys(), [t.phys for t in tets])
        np.testing.assert_array_equal(
            mesh.tri_nodes(),
            [[node_index[c] for c in t.conn[:3]] for t in tris])
        np.testing.assert_array_equal(mesh.tri_phys(), [t.phys for t in tris])

    def test_arrays_are_copies(self, design):
        xyz = design.mesh.node_coords()
        xyz[:] = 0.0
        assert np.any(design.mesh.node_coords() != 0.0)


def _reference_near_field(design, freq, points):
    """The former near_field() search (first tet in mesh order whose
    barycentric coordinates are all >= -1e-10), with each mesh vector read
    once instead of once per vertex."""
    solution = design._last_solutions[freq]
    nodes, tets, node_index = (design.mesh.nodes, design.mesh.tets,
                               design.mesh.nodeIndex)
    E_out = np.zeros((len(points), 3), dtype=complex)
    for pi, point in enumerate(points):
        for tet in tets:
            verts = [nodes[node_index[tet.conn[i]]].xyz for i in range(4)]
            lam = em.compute_barycentric(verts, point)
            if all(v >= -1e-10 for v in lam):
                dofs = np.array([solution[tet.edges[e]] for e in range(6)])
                E_out[pi] = em.evaluate_edge_field(
                    verts, list(tet.edge_orient), dofs, point)
                break
    return E_out


@gmsh_required
def test_near_field_matches_former_search(design):
    design.simulate(2.4e9)
    # The mesh spans [0, 90 mm]^2 in x, y; patch and probe are centred at
    # 45 mm.
    points = np.array([
        [45e-3, 50e-3, 0.8e-3],  # substrate under the patch, off the probe
        [48e-3, 43e-3, 5e-3],    # air above the patch
        [1.0, 1.0, 1.0],         # outside the domain: zero
    ])
    E = design.near_field(2.4e9, points)
    ref = _reference_near_field(design, 2.4e9, points)
    np.testing.assert_allclose(E, ref, rtol=1e-12, atol=0.0)
    assert np.all(np.abs(E[:2]).sum(axis=1) > 0.0)
    assert np.all(E[2] == 0.0)


@gmsh_required
def test_failed_sweep_point_gives_nan_z_and_a_warning(design, monkeypatch):
    """A failed solve (NaN S from the C++ sweep) must not turn into a
    number: Z_in is NaN there, not the 1e6 open-circuit placeholder, and a
    RuntimeWarning names the frequency."""
    real_sweep = em.frequency_sweep
    nan = complex(np.nan, np.nan)

    def sweep_with_failed_middle_point(*args, **kwargs):
        res = real_sweep(*args, **kwargs)
        S = [np.array(m, dtype=complex) for m in res.S_matrices]
        S[1][:, :] = nan

        class Result:
            S_matrices = S
            fields = res.fields
            min_reflection_index = res.min_reflection_index
        return Result()

    monkeypatch.setattr(em, "frequency_sweep", sweep_with_failed_middle_point)
    with pytest.warns(RuntimeWarning, match=r"1 of 3 frequencies.*2\.4 GHz"):
        freqs, z, s11 = design.frequency_sweep(F_START, F_STOP, N_POINTS)
    assert np.isnan(s11[1]) and np.isnan(z[1].real) and np.isnan(z[1].imag)
    assert np.all(np.isfinite(z[[0, 2]]))
    assert np.all(np.abs(z[[0, 2]]) < 1e5)
