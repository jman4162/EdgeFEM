"""VectorEM CAD integration package.

This package provides a light‑weight, modular pipeline for importing
CAD geometry, tagging metadata, generating meshes, previewing results
and exporting a ready‑to‑simulate VectorEM case.  Only a minimal
stub implementation is provided here; the public API mirrors the
intended design so downstream tooling can be developed incrementally."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


@dataclass
class Case:
    """Container holding intermediate artefacts for a CAD case.

    Parameters
    ----------
    geometry : Path | None
        Path to the original CAD geometry.
    mesh : Path | None
        Path to a generated mesh.
    meta : dict
        Optional metadata such as material or boundary tags.
    """

    geometry: Optional[Path] = None
    mesh: Optional[Path] = None
    meta: Dict[str, object] = field(default_factory=dict)


# Convenience imports for SDK style usage.
from .io import step_iges, stl_obj  # noqa: F401
from .mesh.generate import generate_mesh  # noqa: F401
from .ops.tagging import apply_tags  # noqa: F401
from .export.vectorem_case import export_case  # noqa: F401
from .ops.repair import heal_geometry  # noqa: F401

__all__ = [
    "Case",
    "apply_tags",
    "export_case",
    "generate_mesh",
    "heal_geometry",
    "step_iges",
    "stl_obj",
]
