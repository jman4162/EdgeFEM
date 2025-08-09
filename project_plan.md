# VectorEM — Project Plan, README, and Agents

This document contains three sections:

1. **Sprint Plan & Detailed Task Breakdown**
2. **README.md** (ready to drop into the repo root)
3. **AGENTS.md** (design + CI wiring for repo automation agents)

---

## System Design Spec — Platform Update (macOS Apple Silicon)

This amends the earlier System Design Specification:

### 18) Configuration & Deployment (updated)

* **Languages:** C++20, CUDA/HIP (Linux GPU), Python 3.11, Qt 6
* **Platforms:**

  * **Linux (primary)** — full support incl. NVIDIA GPU
  * **Windows (secondary)** — CI builds and tests
  * **macOS Apple Silicon (Tier‑1 local dev)** — verified on **M3 Max**; **CPU-only in v1**. Investigate Metal/MPS backend in v1.x without changing public APIs.
* **Packaging:** macOS local developer setup via Homebrew; `tools/scripts/setup_macos.sh` provided.

### 16) Performance Targets (addendum)

* On an **M3 Max (CPU-only v1)**, the scalar demo completes in **< 60 s**; a small Maxwell case (< \~0.5M dof) completes in **< 5 min** in Release with Ninja. These are smoke-test targets, not guarantees.

### 21) Roadmap (note)

* **v1.x:** Evaluate Metal/MPS-backed iterative kernels (SpMV + smoothers); maintain numerical parity with CPU path.

# 1) Sprint Plan & Detailed Task Breakdown

**Project:** VectorEM (HFSS‑like FEM EM simulator)
**Cadence:** 2‑week sprints (adjust as needed)
**Team Roles:**

* **Numerics Lead** (Maxwell FEM kernels, PML, ports)
* **Mesh/Geometry Lead** (I/O, quality, curvature, boundary layers)
* **Runtime Lead** (solvers, preconditioners, HPC/GPU)
* **DevEx/Infra** (build/CI, packaging, Python SDK, agents)
* **Validation Owner** (benchmarks, MMS, regression packs)

**Global Definition of Done (DoD):**

* Builds on Linux (primary) and Windows (secondary) in CI; **runs locally on macOS Apple Silicon (M1/M2/M3 — verified on M3 Max)** with documented setup
* Unit + integration tests pass; coverage >= 70% on core modules
* Benchmarks match tolerances in Validation Matrix
* User‑level docs/examples updated; CHANGELOG entry present

---

## Sprint 0 — Repo Bootstrap & Scalar Prototype (Weeks 1–2)

**Goals:** Working repo, CI, scalar Helmholtz baseline to prove assembly/solve flow.

**Deliverables:**

* CMake project; Eigen wired; `vectorem_scalar_demo` runs
* Gmsh v2 mesh loader; Tet4 linear element support
* Scalar Helmholtz assembly, basic Dirichlet BCs
* CI (GitHub Actions) with build + unit tests + clang‑tidy

**Tasks:**

* [ ] Init repo structure (`include/`, `src/`, `examples/`, `cmake/`)
* [ ] Add top‑level `CMakeLists.txt` and Eigen finder
* [ ] Implement Gmsh v2 loader (nodes, Tri3, Tet4) + tests
* [ ] Implement scalar Tet4 gradients/volume utilities
* [ ] Assemble scalar Helmholtz (stiffness/mass) + Dirichlet elimination
* [ ] Minimal iterative solver (BiCGSTAB + ILUT) wrapper
* [ ] Example mesh & run script; smoke test prints
* [ ] CI workflow (Linux + Windows), cache dependencies
* [ ] Static analysis: clang‑tidy config and baseline cleanup
* [ ] macOS Apple Silicon (M3 Max) local build smoke test (Xcode clang + Homebrew deps)
* [ ] Add `tools/scripts/setup_macos.sh` and README Mac instructions

**Acceptance:** Demo run prints iterations & 5 solution values; CI green on 2 OSes.

---

## Sprint 1 — Core Maxwell Formulation (Weeks 3–4)

**Goals:** First‑order Nédélec (edge) elements; curl‑curl formulation; complex matrices.

**Deliverables:**

* Edge DOF data structures; local element matrices Ke, Me
* Assembly pipeline `(Ke - ω²Me)e = b` (no ports/PML yet)
* PEC/PMC/symmetry boundary support (prototype)

**Tasks:**

* [ ] Edge (Whitney 1‑form) basis on Tet4; edge indexing map
* [ ] Compute curl(N\_i) and N\_i·N\_j integrals; numerical quadrature
* [ ] Material constants (ε, μ, tanδ) homogeneous; complex arithmetic
* [ ] Boundary ops: PEC (Dirichlet on tangential E), PMC
* [ ] Unit tests: patch vs manufactured solution; eigen cavity sanity

**Acceptance:** Convergence order verified on manufactured problem; cavity eigenfreq within 1%.

---

## Sprint 2 — Ports & Driven‑Modal/Terminal (Weeks 5–6)

**Goals:** Wave/lumped ports; 2D port eigenproblem; S‑parameter extraction path.

**Deliverables:**

* Port cross‑section mesher (Tri/Quad 2D); modal solver (TE/TM modes)
* Port normalization and excitation vector assembly
* S‑parameter computation and Touchstone export

**Tasks:**

* [ ] 2D eigen‑solver (shift‑invert) with PEC/PMC on port boundary
* [ ] Compute modal fields, Z0, power normalization
* [ ] Coupling from port mode to 3D boundary; excitation RHS
* [ ] Lumped ports (terminal) + de‑embedding plane support
* [ ] Post: S‑matrix assembly from fields/port powers; write .sNp
* [ ] Tests: WR‑90 waveguide S11 near cutoff, coax line Z0

**Acceptance:** WR‑90 single‑mode passband shows expected S11/S21 vs analytic; .s2p loads in QUCS/ADS.

---

## Sprint 3 — Radiation & PML (Weeks 7–8)

**Goals:** Open region with robust PML and basic ABC fallback.

**Deliverables:**

* Stretched‑coordinate PML regions with automatic thickness & grading
* PML element Jacobians integrated into Ke/Me
* ABC (1st‑order) optional for coarse runs

**Tasks:**

* [ ] PML region tagging + auto placement from open faces
* [ ] Coordinate stretch functions; stability safeguards
* [ ] Integrate PML into element kernels; unit tests (plane wave absorption)
* [ ] ABC boundary as fallback switch
* [ ] Diagnostics: PML reflection estimates plot

**Acceptance:** Plane wave in box w/ PML shows < −40 dB reflection; patch antenna TRP balance within 1%.

---

## Sprint 4 — Frequency Sweeps & MOR (Weeks 9–10)

**Goals:** Robust sweeps; reuse Krylov/shift data; vector fitting for smooth S(f).

**Deliverables:**

* Discrete per‑frequency driver + caching
* Arnoldi/Ritz reuse between frequencies; optional vector fitting on S(f)
* Sweep policies: robust/balanced/fast

**Tasks:**

* [ ] Frequency driver API + job graph
* [ ] Shift strategies; seed reuse; stopping rules
* [ ] Vector fitting module (stable pole placement); export rational model
* [ ] Tests on filters (monotonic |S21| trends); runtime vs discrete baseline

**Acceptance:** 5× speedup on 401‑pt sweep of 10‑pole filter vs naive discrete solves.

---

## Sprint 5 — Adaptive h/p & Error Estimators (Weeks 11–12)

**Goals:** Goal‑oriented refinement toward S‑param/far‑field functionals.

**Deliverables:**

* Residual‑based and goal‑oriented estimators; element marking
* h‑refine (1‑irregular) and p‑raise (orders 1→3); hp policy heuristic

**Tasks:**

* [ ] Residual estimator for curl‑curl; data structures for element error
* [ ] Dual solve for goal functional (S\_ij or TRP) (coarse adjoint)
* [ ] Mesh refinement (local tetra split) + conformity fixes
* [ ] p‑enrichment based on smoothness indicator
* [ ] Stop criteria: |ΔS|, estimator threshold, wallclock cap

**Acceptance:** Patch antenna gain converges to within 0.5 dB with < 60% dof vs uniform refine.

---

## Sprint 6 — Post‑Processing: Fields & Patterns (Weeks 13–14)

**Goals:** Near‑field visualization; near‑to‑far; antenna metrics.

**Deliverables:**

* Field probes, cut‑planes, Poynting vector
* Huygens surface extraction; far‑field patterns, gain, AR, cross‑pol
* Pattern exports (.pat JSON/CSV); Smith/Bode plotting

**Tasks:**

* [ ] Field sampling/interpolation on cut planes and probes
* [ ] Huygens equivalent surface and Stratton–Chu variant
* [ ] Polar/3D plots; integration for TRP/efficiency
* [ ] Report generator templates

**Acceptance:** Patch antenna 2.45 GHz pattern within literature values (beamwidth, F/B, efficiency).

---

## Sprint 7 — Runtime & Solvers (Weeks 15–16)

**Goals:** Robust linear algebra; domain decomposition; preconditioning.

**Deliverables:**

* GMRES/FGMRES + auxiliary‑space AMG for Maxwell
* Multifrontal direct solver integration (optional vendor API)
* Domain decomposition preconditioner (BDD/FETI‑like, coarse space)

**Tasks:**

* [ ] Auxiliary space preconditioner implementation and benchmarks
* [ ] Direct solver wrapper + pivoting options
* [ ] Stress tests on ill‑conditioned geometries

**Acceptance:** 10M dof case solves with < 30 GMRES iterations using AMG precond (lab benchmark).

---

## Sprint 8 — HPC, GPU, and Checkpointing (Weeks 17–18)

**Goals:** Hybrid MPI+threads; optional GPU SpMV/AMG; fault‑tolerant restarts.

**Deliverables:**

* MPI parallel assembly/solve path; SLURM integration
* Optional CUDA/HIP backends for SpMV and smoothers
* Checkpointing for factorization/Krylov basis

**Tasks:**

* [ ] Partitioning (ParMETIS) + ghost exchange for assembly
* [ ] MPI collectives & overlap; mixed precision experiments
* [ ] GPU kernels for hot paths; fallback when unavailable
* [ ] Checkpoint/restart files; versioning

**Acceptance:** Strong scaling ≥ 70% to 64 ranks; 2–4× GPU speedup on iterative cases.

---

## Sprint 9 — Python SDK & CLI Polishing (Weeks 19–20)

**Goals:** Scriptability; headless runs; notebooks for validation.

**Deliverables:**

* pybind11 module (`pyvectorem`); high‑level Python API
* CLI subcommands (mesh/solve/post)
* Jupyter notebooks: S‑params, patterns, convergence studies

**Tasks:**

* [ ] Stable Python API surface and docstrings
* [ ] Touchstone & HDF5 readers/writers
* [ ] Examples repository and CI notebook smoke tests

**Acceptance:** Two end‑to‑end examples runnable via Python and CLI.

---

## Sprint 10 — UX/GUI (Weeks 21–22)

**Goals:** Minimal Qt GUI mimicking HFSS tree and inspectors.

**Deliverables:**

* Project tree (Geometry/Setup/Excitations/Analysis/Results)
* Property panels; live plot widgets; wizards for TL/patch/waveguide

**Tasks:**

* [ ] Qt 6 app skeleton; async job runner
* [ ] Mesh/port preview; PML grading plot
* [ ] Report viewer and export to PNG/PDF

**Acceptance:** Demo video: create waveguide model, define port, sweep, view S11 plot.

---

## Sprint 11 — Validation, QA, and Release (Weeks 23–24)

**Goals:** Golden benchmarks; docs; packaging.

**Deliverables:**

* Validation pack (20 cases) + tolerances; nightly regression
* User docs site; CHANGELOG; LICENSE
* Installers (Linux) and wheels for Python SDK

**Tasks:**

* [ ] MMS tests; cavity, TLs, WG steps, patch antenna, filter
* [ ] Docs site (MkDocs or Sphinx) with tutorials
* [ ] Packaging scripts; versioning; release checklist

**Acceptance:** v1.0 tag cut; artifacts uploaded; docs published; demo webinar deck.

---

## Cross‑Cutting Backlog (Prioritize as needed)

* [ ] Conductor roughness models (Huray/Groiss); temperature coeffs
* [ ] Stackup/PCB import (ODB++, Gerber + layer stack)
* [ ] Debye/Lorentz dispersive materials; passivity enforcement
* [ ] Improved MOR; passivity‑preserving rational fits
* [ ] Pattern mask compliance checks (regulatory)
* [ ] Plugin system for custom post‑processors

## Risks & Mitigations

* **PML instability:** auto tuning + validation cases; ABC fallback
* **AMG robustness:** auxiliary‑space design + conservative defaults; direct solver fallback
* **Mesh quality:** curved elements + defeaturing; sliver detection & repair
* **User port setup errors:** port field preview + orthogonality checks + wizards

---

