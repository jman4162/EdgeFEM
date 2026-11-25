"""Binding-level input handling: numpy strides, sizes, dtypes."""

import numpy as np
import pytest

em = pytest.importorskip("pyedgefem", reason="pyedgefem not built")


def test_vecc_from_strided_array():
    """A strided view (arr[::2]) must give its own elements, not the
    underlying buffer's first n."""
    base = np.arange(10) * (1 + 1j)
    v = em.VecC(base[::2])
    assert len(v) == 5
    assert [v[i] for i in range(5)] == list(base[::2])


def test_vecc_from_real_array():
    v = em.VecC(np.array([1.0, 2.0, 3.0]))
    assert [v[i] for i in range(3)] == [1, 2, 3]


def _unit_tet():
    return [np.array(p, dtype=float) for p in
            ([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])]


def test_evaluate_edge_field_checks_dof_count():
    with pytest.raises(ValueError, match="6"):
        em.evaluate_edge_field(_unit_tet(), [1] * 6,
                               np.ones(5, dtype=complex),
                               np.array([0.25, 0.25, 0.25]))


def test_evaluate_edge_field_strided_dofs():
    dofs = np.arange(12, dtype=complex)
    a = em.evaluate_edge_field(_unit_tet(), [1] * 6, dofs[::2],
                               np.array([0.2, 0.3, 0.1]))
    b = em.evaluate_edge_field(_unit_tet(), [1] * 6, dofs[::2].copy(),
                               np.array([0.2, 0.3, 0.1]))
    np.testing.assert_allclose(a, b)


def test_port_abc_scale_docstring_is_not_stale():
    doc = em.MaxwellParams.port_abc_scale.__doc__ or ""
    assert "0.5 optimal" not in doc


def _wr90():
    mesh = em.load_gmsh("examples/rect_waveguide.msh")
    bc = em.build_edge_pec(mesh, 1)
    dims = em.RectWaveguidePort()
    dims.a, dims.b = 0.02286, 0.01016
    mode = em.solve_te10_mode(dims, 10e9)
    kc_sq = (np.pi / dims.a) ** 2
    ports = [em.build_wave_port_2d(mesh, t, mode, set(bc.dirichlet_edges),
                                   kc_sq) for t in (2, 3)]
    p = em.MaxwellParams()
    p.omega = 2 * np.pi * 10e9
    return mesh, p, bc, ports


def test_solves_release_the_gil():
    """While a solve runs in another thread, the main thread keeps running
    Python: the longest pause it sees is well below one solve."""
    import threading
    import time
    mesh, p, bc, ports = _wr90()
    t0 = time.perf_counter()
    em.calculate_sparams_eigenmode(mesh, p, bc, ports)
    one_solve = time.perf_counter() - t0

    done = threading.Event()

    def work():
        for _ in range(3):
            em.calculate_sparams_eigenmode(mesh, p, bc, ports)
        done.set()

    th = threading.Thread(target=work)
    stamps = []
    th.start()
    while not done.is_set():
        stamps.append(time.perf_counter())
    th.join()
    gaps = np.diff(stamps)
    assert len(gaps) > 0
    assert gaps.max() < 0.5 * one_solve, (gaps.max(), one_solve)


def test_embedded_patterns_are_not_implemented():
    """The embedded-pattern functions summed raw edge DOFs of the whole
    mesh (no Huygens surface); they refuse to run until a real
    implementation exists."""
    mesh, p, bc, ports = _wr90()
    cfg = em.ArrayPatternConfig()
    cfg.theta_rad = [0.0, 0.5]
    cfg.k0 = 2 * np.pi * 10e9 / 299792458.0
    with pytest.raises(NotImplementedError, match="extract_huygens_surface"):
        em.compute_embedded_patterns(mesh, p, bc, ports, cfg)
    with pytest.raises(NotImplementedError, match="extract_huygens_surface"):
        em.compute_single_embedded_pattern(mesh, p, bc, ports, 0, cfg)


def test_package_all_is_consistent():
    """Every name in edgefem.__all__ exists, the documented entry points are
    exported, and the deprecated normalize_port_weights is not."""
    import edgefem
    missing = [n for n in edgefem.__all__ if not hasattr(edgefem, n)]
    assert not missing, missing
    for name in ("PortKind", "build_lumped_port", "LumpedPortConfig",
                 "build_wave_port_2d", "frequency_sweep",
                 "calculate_sparams_with_fields",
                 "calculate_sparams_eigenmode_with_fields",
                 "extract_huygens_surface"):
        assert name in edgefem.__all__, name
    assert "normalize_port_weights" not in edgefem.__all__
    assert "110 GHz" not in (edgefem.__doc__ or "")
