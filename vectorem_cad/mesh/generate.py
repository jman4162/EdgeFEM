"""Mesh generation utilities using Gmsh."""

from __future__ import annotations

from pathlib import Path

from .. import Case


def generate_mesh(case: Case, element_size: float = 1.0) -> Case:
    """Generate a tetrahedral mesh using Gmsh.

    This stub merely records a fake mesh file path derived from the
    geometry.  Real implementations would invoke the Gmsh API and write
    a ``.msh`` file containing physical groups.
    """

    if not case.geometry:
        raise ValueError("Case has no geometry to mesh")

    mesh_path = Path(str(case.geometry) + ".msh")
    case.mesh = mesh_path
    return case


__all__ = ["generate_mesh"]
