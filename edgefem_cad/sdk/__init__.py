"""Public SDK convenience layer.

This module re-exports the high level functions so that users can write
``from edgefem_cad.sdk import import_cad`` etc.  The goal is to keep the
public surface tiny and stable."""

from .. import Case, apply_tags, export_case, generate_mesh
from ..io.step_iges import load_geometry as import_cad

__all__ = [
    "Case",
    "apply_tags",
    "export_case",
    "generate_mesh",
    "import_cad",
]
