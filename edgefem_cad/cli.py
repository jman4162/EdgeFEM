"""Command line interface for :mod:`edgefem_cad`.

The CLI is intentionally minimal yet mirrors the proposed commands in
`system_design_spec.md`.  Each sub‑command delegates to the stub API
functions defined in this package.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from . import Case
from .export.edgefem_case import export_case
from .io import step_iges, stl_obj
from .mesh.generate import generate_mesh
from .ops.repair import heal_geometry
from .ops.tagging import apply_tags
from .viz.preview import preview


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="edgefem-cad")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_import = sub.add_parser("import", help="Import a CAD file")
    p_import.add_argument("file", help="Input CAD file")
    p_import.add_argument("-o", "--out", default="case", help="Output directory")

    p_mesh = sub.add_parser("mesh", help="Generate mesh")
    p_mesh.add_argument("case_dir")

    p_tag = sub.add_parser("tag", help="Apply tagging")
    p_tag.add_argument("case_dir")

    p_preview = sub.add_parser("preview", help="Preview geometry or mesh")
    p_preview.add_argument("case_dir")
    p_preview.add_argument("--mesh", action="store_true")

    p_export = sub.add_parser("export", help="Export EdgeFEM case")
    p_export.add_argument("case_dir")
    p_export.add_argument("--format", default="edgefem")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "import":
        path = Path(args.file)
        if path.suffix.lower() in step_iges.SUPPORTED_SUFFIXES:
            case = step_iges.load_geometry(path)
        else:
            case = stl_obj.load_mesh(path)
        export_case(case, args.out)
        return 0

    if args.cmd == "mesh":
        case = Case(geometry=Path(args.case_dir) / "geometry.step")
        generate_mesh(case)
        return 0

    if args.cmd == "tag":
        case = Case()
        apply_tags(case)
        return 0

    if args.cmd == "preview":
        case = Case()
        preview(case, show_mesh=args.mesh)
        return 0

    if args.cmd == "export":
        case = Case(mesh=Path(args.case_dir) / "mesh.msh")
        export_case(case, args.case_dir)
        return 0

    parser.error(f"Unknown command {args.cmd}")
    return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
