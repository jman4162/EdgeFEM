# AGENTS.md — Developer Tooling & CI Automation

This repo includes small Python CLIs under `tools/` that automate routine
development tasks. They run locally and in CI (`.github/workflows/ci.yml`).
None of them auto-commit; they produce artifacts for maintainers to inspect.

## What exists today

### Mesh QA (`tools/agents/mesh_qa.py`)
Parses Gmsh v2 ASCII meshes and reports per-mesh quality metrics as JSON:
node/element counts, tetrahedron signed volumes (degenerate/inverted
detection), and edge-length aspect ratios. Mirrors the checks in
`src/mesh_quality.cpp`.

```bash
python tools/agents/mesh_qa.py --in examples --out artifacts/mesh_qa
```

CI runs this over `examples/` on every push/PR and uploads the JSON report.

### Simulation scaffolder (`tools/new_simulation.py`)
Generates a `.geo` geometry file and a `run_simulation.py` driver from a
template (waveguide, dipole, unit cell) into `simulations/<name>/`.

### Pattern plotter (`tools/plot_pattern.py`)
Polar radiation-pattern plots from CSV pattern exports.

## CI

The real pipeline is `.github/workflows/ci.yml`: build (Ninja, Release,
`EDGEFEM_PYTHON=ON`), clang-format/clang-tidy via `tools/lint.sh`, `ctest`,
`pytest python/tests/`, then mesh QA. `wheels.yml` builds wheels on
Linux/macOS; `docs.yml` deploys MkDocs.

## Ideas not yet implemented

Earlier drafts of this file described a larger agent catalog (numerics PR
review, docs sync, benchmark triage, release-notes generation) and an MCP
server for repository retrieval. None of that exists yet. If you add such a
tool, follow the conventions below rather than reintroducing aspirational
docs.

## Adding a new tool

1. Create `tools/agents/<name>.py` with `--in`/`--out` args and JSON output.
2. Add a CI step in `.github/workflows/ci.yml` and mention it here.
3. Include a sample input and expected output in `tools/agents/samples/`.
