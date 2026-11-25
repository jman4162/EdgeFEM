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
