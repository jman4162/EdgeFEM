# VectorEM — Codex Tasks

## How to use
At the start of each Codex session:
1. Open and read `README.md`, `AGENTS.md`, and `tools/scripts/setup_macos.sh`.
2. Follow repo rules, build, and run smoke tests:
   ```bash
   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ctest --test-dir build -L smoke -j

3. Return a **unified diff**, **commit message**, and the **full terminal log**.

---

## Meta: Repo rules & workflow (run first)

**Prompt:**

> You are contributing to the **VectorEM** repo. Before doing anything, read `README.md`, `AGENTS.md`, and `tools/scripts/setup_macos.sh`. Do not change public APIs unless tests and docs are updated in the same PR. Prefer small diffs. After edits, run the exact build/test commands above. Return a unified diff, a concise commit message, and the full configure/build/test logs.

---

## Sprint‑0 / Infra

### T1 — Bootstrap CI + smoke pipeline

**Goal:** Add GitHub Actions CI (`.github/workflows/ci.yml`) that builds, runs `ctest -L smoke`, runs a stub mesh QA agent, and uploads artifacts.
**Edits:**

* Add workflow as specified in **AGENTS.md**.
* Create a no‑op `tools/agents/mesh_qa.py` that prints a small JSON summary.
  **DoD:** CI is green on Ubuntu; local `ctest -L smoke` passes.

### T2 — Example mesh so smoke test is turnkey

**Goal:** Provide a tiny cavity mesh.
**Edits:**

* Add `examples/cube_cavity.geo` and script `examples/make_mesh.sh` to generate `examples/cube_cavity.msh` (Gmsh v2 tetra). Ensure exterior Tri faces carry `phys=1`.
* Add `examples/README.md` with generation steps.
  **DoD:** `vectorem_scalar_demo examples/cube_cavity.msh 1` finishes < 60 s on macOS; `ctest -L smoke` passes.

### T3 — clang‑tidy + format guardrails

**Goal:** Enforce style and static analysis.
**Edits:**

* Add `.clang-tidy` and `.clang-format` (LLVM base, 2‑space indent).
* Add `tools/lint.sh` and a CI step named **Static analysis** that fails on new warnings.
  **DoD:** CI fails if a PR adds new tidy warnings in `src/`.

---

## Numerics Core

### T4 — Edge DOF scaffolding (Nédélec Tet(1))

**Goal:** Introduce edge DOFs and local basis utilities.
**Edits:**

* Add `include/vectorem/edge_basis.hpp` and `src/edge_basis.cpp` with Whitney 1‑forms on Tet(1), routines for `curl(N_i)` and integrals `∫ curl(N_i)·curl(N_j) dV` and `∫ N_i·N_j dV`.
* Add a global edge indexing map and element→edge connectivity.
* Keep scalar prototype behind a flag; both paths must build.
  **DoD:** New unit test `tests/test_edge_indexing.cpp` validates consistent edge numbering and local→global mappings; smoke still passes.

### T5 — Maxwell curl–curl assembly skeleton

**Goal:** Implement `(Ke − ω²Me)e = b` for homogeneous ε, μ.
**Edits:**

* Add `src/assemble_maxwell.cpp` (+ header) to assemble complex sparse A and RHS.
* Add PEC boundary handling for tangential E on tagged Tri faces.
  **DoD:** Manufactured‑solution or cavity sanity unit test compiles and runs < 60 s; scalar path unaffected.

---

## Ports & S‑parameters

### T6 — 2D port eigenmode (prototype)

**Goal:** Solve port cross‑section eigenmodes.
**Edits:**

* Add `src/ports/port_eigensolve.cpp` and headers in `include/vectorem/ports/`.
* Minimal 2D triangular mesher or load a 2D mesh; shift‑invert eigen solve; compute `Z0` and 1 W normalization.
  **DoD:** WR‑90 port example outputs cutoff freq and `Z0` close to analytic.

### T7 — Couple port mode to 3D and compute S

**Goal:** Boundary excitation and S‑matrix.
**Edits:**

* Implement coupling of port mode to 3D boundary to form RHS.
* Add Touchstone writer `src/io/touchstone.cpp`.
  **DoD:** Straight WR‑90 section shows in‑band `|S11| < −20 dB`; writes `out/test.s2p`.

---

## Open Boundaries

### T8 — PML (stretched coordinates)

**Goal:** Open region modeling.
**Edits:**

* Tag PML regions; implement coordinate stretch; modify Ke/Me assembly inside PML.
* Add ABC (1st‑order) fallback flag.
  **DoD:** Plane‑wave absorption test ≤ −40 dB reflection; add `ctest -L pml` (skipped by default in CI).

---

## Sweeps & Speed

### T9 — Frequency sweep driver + reuse

**Goal:** Faster sweeps.
**Edits:**

* Implement a sweep controller supporting discrete, seed reuse, and a “balanced” policy.
* CLI option to choose sweep policy; timing logs.
  **DoD:** On a toy 2‑resonator example, 401‑pt sweep is ≥ 2× faster than naive discrete.

### T10 — Vector fitting of S(f)

**Goal:** Smooth models.
**Edits:**

* Add `src/mor/vector_fit.cpp` to fit S(f) with passivity checks; export a JSON model.
  **DoD:** Fit RMS error < 1% on the toy sweep; model file written.

---

## Post‑Processing

### T11 — Near‑to‑far (Huygens surface)

**Goal:** Far‑field patterns.
**Edits:**

* Implement equivalent current extraction on a Huygens surface and Stratton–Chu far‑field integration in `src/post/ntf.cpp`.
* Add `tools/plot_pattern.py` to visualize CSV patterns.
  **DoD:** Patch antenna example produces 2D polar pattern and CSV export.

---

## Python / DevEx

### T12 — pybind11 bindings v1

**Goal:** Scriptability.
**Edits:**

* Build `pyvectorem` exposing: mesh load, scalar assemble/solve, and Touchstone export.
* Add a tiny Python smoke test runnable in CI (< 30 s).
  **DoD:** `python -c "import pyvectorem as em; print('ok')"` works; CI job runs the Python smoke test.

---

## Prompt Template (reuse)

> **Context:** VectorEM repository. Read `README.md`, `AGENTS.md`, and `tools/scripts/setup_macos.sh` first.
> **Task:** \<one clear, small objective>
> **Constraints:** Keep public APIs stable; minimal diffs; add/adjust tests; support macOS Apple Silicon (CPU) and Linux.
> **Commands:**
>
> ```bash
> cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
> cmake --build build -j
> ctest --test-dir build -L smoke -j
> ```
>
> **Deliverables:** Unified diff, commit message, full terminal logs; include full contents for any new files.
> **Definition of Done:** \<explicit pass/fail checks, perf/time bounds, file artifacts>

```
```
