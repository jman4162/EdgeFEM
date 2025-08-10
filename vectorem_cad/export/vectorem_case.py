"""Export a :class:`~vectorem_cad.Case` into a VectorEM case directory."""

from __future__ import annotations

import json
from pathlib import Path

from .. import Case


def export_case(case: Case, out_dir: str | Path) -> Path:
    """Write a minimal VectorEM case structure.

    The function creates ``mesh.msh`` and ``settings.json`` placeholders
    so that downstream tooling has consistent expectations.
    """

    dest = Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)

    if case.mesh:
        # In the stub we do not actually copy any mesh data.
        (dest / "mesh.msh").write_text("0 $MeshFormat\n4.1 0 8\n0\n$EndMeshFormat\n")

    settings = {"notes": "Placeholder case"}
    (dest / "settings.json").write_text(json.dumps(settings, indent=2))
    return dest


__all__ = ["export_case"]
