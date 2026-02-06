"""Visualization helpers."""

from __future__ import annotations

from pathlib import Path

from .. import Case


def preview(case: Case, show_mesh: bool = False) -> Path:
    """Return a path to a dummy screenshot for tests.

    The actual project would use :mod:`pyvista` or ``gmsh`` for
    interactive viewing.  Here we just create an empty PNG file to
    simulate an exported screenshot.
    """

    screenshot = Path("preview.png")
    screenshot.write_bytes(b"")
    return screenshot


__all__ = ["preview"]
