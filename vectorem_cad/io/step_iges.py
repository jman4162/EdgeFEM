"""Import STEP and IGES geometry using `pythonocc-core`.

The real implementation would convert BREP data into an internal
representation.  For now we only record the file path and perform
basic suffix validation."""

from pathlib import Path
from .. import Case

SUPPORTED_SUFFIXES = {".step", ".stp", ".iges", ".igs"}


def load_geometry(path: str | Path) -> Case:
    """Load a STEP or IGES file and return a :class:`Case`.

    Parameters
    ----------
    path:
        Location of the CAD file.

    Returns
    -------
    Case
        The newly created case with the geometry path set.
    """

    file_path = Path(path)
    if file_path.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise ValueError(f"Unsupported CAD format: {file_path.suffix}")
    return Case(geometry=file_path)


__all__ = ["load_geometry", "SUPPORTED_SUFFIXES"]
