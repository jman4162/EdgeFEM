#!/usr/bin/env python3
"""Mesh QA agent.

Parses Gmsh v2 ASCII meshes and reports quality metrics as JSON:
node/element counts, tetrahedron signed volumes (degenerate/inverted
detection), and edge-length aspect ratios. Mirrors the checks in
src/mesh_quality.cpp.
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

DEGENERATE_VOLUME = 1e-15
ASPECT_WARN = 20.0


def parse_msh_v2(path):
    """Return (nodes, tets) from a Gmsh v2 ASCII file.

    nodes: dict node_id -> (x, y, z); tets: list of 4-tuples of node ids.
    """
    nodes = {}
    tets = []
    with open(path) as f:
        lines = iter(f)
        for line in lines:
            tag = line.strip()
            if tag == "$MeshFormat":
                version = next(lines).split()[0]
                if not version.startswith("2"):
                    raise ValueError(
                        f"Gmsh format {version} unsupported; EdgeFEM requires "
                        "v2 ASCII (gmsh -format msh2)")
            elif tag == "$Nodes":
                count = int(next(lines))
                for _ in range(count):
                    parts = next(lines).split()
                    nodes[int(parts[0])] = tuple(float(v) for v in parts[1:4])
            elif tag == "$Elements":
                count = int(next(lines))
                for _ in range(count):
                    parts = next(lines).split()
                    elem_type = int(parts[1])
                    if elem_type == 4:  # Tet4
                        num_tags = int(parts[2])
                        conn = parts[3 + num_tags:7 + num_tags]
                        tets.append(tuple(int(n) for n in conn))
    return nodes, tets


def tet_signed_volume(p0, p1, p2, p3):
    a = [p1[i] - p0[i] for i in range(3)]
    b = [p2[i] - p0[i] for i in range(3)]
    c = [p3[i] - p0[i] for i in range(3)]
    det = (a[0] * (b[1] * c[2] - b[2] * c[1])
           - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    return det / 6.0


def tet_edge_aspect(points):
    lengths = []
    for i, j in itertools.combinations(range(4), 2):
        d = sum((points[i][k] - points[j][k]) ** 2 for k in range(3)) ** 0.5
        lengths.append(d)
    lo = min(lengths)
    return float("inf") if lo == 0.0 else max(lengths) / lo


def qa_mesh(path):
    nodes, tets = parse_msh_v2(path)
    report = {
        "file": str(path),
        "num_nodes": len(nodes),
        "num_tets": len(tets),
    }
    if not tets:
        report["status"] = "no_tets"
        return report

    volumes = []
    aspects = []
    inverted = degenerate = 0
    for conn in tets:
        pts = [nodes[n] for n in conn]
        vol = tet_signed_volume(*pts)
        volumes.append(vol)
        aspects.append(tet_edge_aspect(pts))
        if abs(vol) < DEGENERATE_VOLUME:
            degenerate += 1
        elif vol < 0:
            inverted += 1

    report.update({
        "total_volume": sum(abs(v) for v in volumes),
        "min_abs_volume": min(abs(v) for v in volumes),
        "max_aspect_ratio": max(aspects),
        "num_inverted": inverted,
        "num_degenerate": degenerate,
        "num_high_aspect": sum(1 for a in aspects if a > ASPECT_WARN),
    })
    report["status"] = "ok" if inverted == 0 and degenerate == 0 else "fail"
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Gmsh v2 mesh QA agent")
    parser.add_argument("--in", dest="in_path", required=True,
                        help="Mesh file or directory containing .msh files")
    parser.add_argument("--out", dest="out_dir", required=True,
                        help="Output directory for summary.json")
    args = parser.parse_args()

    in_path = Path(args.in_path)
    mesh_files = sorted(in_path.glob("*.msh")) if in_path.is_dir() else [in_path]

    reports = []
    for mesh in mesh_files:
        try:
            reports.append(qa_mesh(mesh))
        except Exception as exc:  # malformed mesh should fail QA, not crash CI
            reports.append({"file": str(mesh), "status": "parse_error",
                            "error": str(exc)})

    failed = [r["file"] for r in reports if r["status"] not in ("ok", "no_tets")]
    summary = {"num_meshes": len(reports), "failed": failed, "meshes": reports}

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"num_meshes": len(reports), "failed": failed}))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
