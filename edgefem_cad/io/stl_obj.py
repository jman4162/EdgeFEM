"""Import STL and OBJ meshes using :mod:`trimesh` if available."""

from pathlib import Path
from .. import Case

SUPPORTED_SUFFIXES = {".stl", ".obj"}


def load_mesh(path: str | Path) -> Case:
    """Load a facet mesh and return a :class:`Case` instance.

    This placeholder simply verifies the suffix and records the path.
    """

    file_path = Path(path)
    if file_path.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise ValueError(f"Unsupported mesh format: {file_path.suffix}")
    return Case(geometry=file_path)


__all__ = ["load_mesh", "SUPPORTED_SUFFIXES"]
