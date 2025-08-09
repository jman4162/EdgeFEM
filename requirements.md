# REQUIREMENTS.md — VectorEM v1.0

**Status:** Draft (v0.1)  
**Owner:** Core Team (Numerics Lead + Runtime Lead)  
**Last updated:** <set on commit>

## 1. Purpose & Scope
VectorEM is a 3D frequency‑domain FEM electromagnetics solver for RF/mmWave analysis. v1.0 targets accurate S‑parameters and antenna radiation patterns with first‑order Nédélec elements, wave/lumped ports, and PML. This document defines testable **functional and non‑functional requirements** for v1.0.

## 2. Definitions & Acronyms
- **FEM:** Finite Element Method  
- **Nédélec(1):** First‑order edge elements (Whitney 1‑forms)  
- **PML:** Perfectly Matched Layer (stretched‑coordinate)  
- **PEC/PMC:** Perfect Electric/Perfect Magnetic Conductor  
- **S‑parameters:** Network scattering parameters exported as Touchstone `.sNp`

## 3. Assumptions & Dependencies
- C++20 toolchain; CMake ≥ 3.20; Eigen ≥ 3.4  
- Linux is the primary platform; Windows used for CI validation; macOS Apple Silicon used for local development (CPU‑only in v1).  
- GPU acceleration (CUDA/HIP) is Linux‑only in v1. macOS Metal/MPS research may land in v1.x without API changes.

## 4. Functional Requirements (FR)
Each item is labeled with an ID; “MoSCoW” priority in brackets. **Verification** method shown at end.

### Geometry, Meshing, Materials
- **FR‑001 [Must]:** Load Gmsh v2/v4 tetrahedral meshes with physical tags for regions and boundary faces.  
  _Verify:_ unit/integration tests with sample meshes.
- **FR‑002 [Should]:** Import STEP → internal mesh via external mesher pipeline (deferred if not ready).  
  _Verify:_ demo workflow doc; optional test.
- **FR‑003 [Must]:** Assign material properties per region: ε_r(f), μ_r(f) (constant in v1), σ, tanδ; allow complex ε, μ.  
  _Verify:_ unit tests on element assembly with complex params.
- **FR‑004 [Should]:** Conductor surface impedance boundary (SIBC) with frequency‑dependent σ; roughness models may be added in v1.x.  
  _Verify:_ boundary unit tests.

### Maxwell Solver Core
- **FR‑010 [Must]:** Assemble curl–curl system with Nédélec(1) on Tet4: `(Ke − ω² Me) e = b`, complex sparse.  
  _Verify:_ manufactured‑solution convergence test; cavity eigen sanity.
- **FR‑011 [Must]:** Boundary conditions: PEC, PMC, symmetry plane (tangential E handling).  
  _Verify:_ unit tests on tagged faces.
- **FR‑012 [Must]:** Iterative linear solver (GMRES/FGMRES) with auxiliary‑space AMG preconditioner; fallback to direct solver where available.  
  _Verify:_ solver regression on medium problem; iteration bound.

### Ports & Network Quantities
- **FR‑020 [Must]:** Wave ports with 2D eigenmode solve (fundamental + optional higher modes), power‑normalized to 1 W; compute mode impedances.  
  _Verify:_ WR‑90 example vs analytic cutoff and Z0.
- **FR‑021 [Must]:** Lumped/terminal ports; optional series/parallel RLC.  
  _Verify:_ coax/strip example and de‑embedding test.
- **FR‑022 [Must]:** S‑parameter extraction and export to Touchstone `.sNp` with frequency list.  
  _Verify:_ file schema check; cross‑tool import.

### Open Boundaries & Radiation
- **FR‑030 [Must]:** PML regions (auto placement from open faces + manual tagging) using stretched coordinates.  
  _Verify:_ plane‑wave absorption case (≤ −40 dB reflection).
- **FR‑031 [Should]:** ABC (1st‑order) fallback for coarse passes.  
  _Verify:_ test toggling PML/ABC; numerical stability check.
- **FR‑032 [Must]:** Near‑to‑far transformation via Huygens surface (equivalent currents) to compute far‑field patterns (gain, directivity, AR, cross‑pol).  
  _Verify:_ patch antenna pattern vs literature; TRP balance.

### Sweeps, Adaptivity, and Post
- **FR‑040 [Must]:** Frequency sweep controller supporting discrete solves and seed/Krylov reuse policy (“balanced”).  
  _Verify:_ runtime improvement on toy model vs discrete baseline.
- **FR‑041 [Should]:** Vector fitting of S(f) for smooth models; passivity checks (clamp or report).  
  _Verify:_ RMS fit error threshold on sample sweep.
- **FR‑042 [Should]:** Goal‑oriented h/p adaptivity (element marking by residual estimator; p‑raise policy); can be flagged experimental.  
  _Verify:_ convergence to tolerance on patch/filter with fewer DOFs than uniform refine.
- **FR‑043 [Must]:** Field post‑processing: |E|, |H|, Poynting, cut‑planes, probes; CSV/HDF5 export.  
  _Verify:_ sample outputs and schema tests.

### Interfaces & UX
- **FR‑050 [Must]:** CLI to run: mesh check, solve, post; exit codes on failure.  
  _Verify:_ CLI tests.
- **FR‑051 [Must]:** Python SDK (`pyvectorem`) exposing: mesh load, solve, S‑params export, basic field sampling.  
  _Verify:_ Python smoke test in CI (< 30 s).
- **FR‑052 [Should]:** Minimal Qt GUI (tree view, property panels, plots) for demos.  
  _Verify:_ manual demo script and screenshot artifact.

## 5. Non‑Functional Requirements (NFR)

### Accuracy & Numerical Quality
- **NFR‑100 [Must]:** Step‑to‑step mesh‑refinement delta on S‑params **|ΔS| ≤ 0.01** before stop.  
  _Verify:_ adaptive runs store ΔS logs.
- **NFR‑101 [Must]:** Cavity eigenfrequency error **≤ 0.5%** vs analytic for low‑order modes.  
  _Verify:_ eigen sanity test.
- **NFR‑102 [Must]:** Far‑field power balance (TRP) error **≤ 1%** for closed Huygens surface.  
  _Verify:_ post‑processing check.
- **NFR‑103 [Must]:** PML total reflection **≤ −40 dB** on plane‑wave incidence benchmark.  
  _Verify:_ absorption test.
- **NFR‑104 [Should]:** Port mode orthogonality leakage **≤ −30 dB** between distinct modes.  
  _Verify:_ port solver test.

### Performance & Scalability (targets; not hard fails unless stated)
- **NFR‑120 [Must]:** macOS **M3 Max (CPU‑only)** smoke test completes `< 60 s` on the tiny scalar case; small Maxwell case (`< ~0.5M` DOF) completes `< 5 min` Release build.  
  _Verify:_ `ctest -L smoke`, timed Maxwell demo.
- **NFR‑121 [Should]:** 8‑core workstation solves a **2–5M DOF** frequency point in `< 30 min` (iterative) given a well‑conditioned case and robust preconditioning.  
  _Verify:_ lab benchmark.
- **NFR‑122 [Should]:** On a 64‑rank cluster, strong scaling efficiency **≥ 70%** up to **50M DOF** assemblies.  
  _Verify:_ scaling benchmark.

### Reliability & Operability
- **NFR‑140 [Must]:** Deterministic runs: hash of CAD+mesh+settings recorded; results annotated with `run_id`.  
  _Verify:_ manifest file present; reproducibility spot check.
- **NFR‑141 [Must]:** Checkpoint/restart for long jobs (factorization/Krylov basis) with versioned files.  
  _Verify:_ restart test on interrupted run.
- **NFR‑142 [Must]:** Structured logging with levels; default INFO, solver timings at DEBUG.  
  _Verify:_ logs contain stage timings.

### Security & Privacy
- **NFR‑160 [Must]:** Solvers run with **no external network access** by default (HPC‑safe).  
  _Verify:_ CI test mocks network; runtime flag.
- **NFR‑161 [Should]:** Optional telemetry is **opt‑in** and anonymized; off by default.  
  _Verify:_ config switch; data schema doc.
- **NFR‑162 [Should]:** Plugins are checksummed/signed; untrusted plugins blocked by default.  
  _Verify:_ loader test.

### Platform & Compatibility
- **NFR‑180 [Must]:** Build & test in CI on **Linux (primary)** and **Windows (secondary)**.  
  _Verify:_ CI green on both.
- **NFR‑181 [Must]:** **Run locally on macOS Apple Silicon (M1/M2/M3)** with documented setup; **verified on M3 Max**.  
  _Verify:_ `tools/scripts/setup_macos.sh` and smoke test pass.
- **NFR‑182 [Should]:** Export **Touchstone v1.1** compliant `.sNp` files.  
  _Verify:_ schema check; import into third‑party tool.

### Code Quality & Docs
- **NFR‑200 [Must]:** Unit + integration test coverage **≥ 70%** on core modules.  
  _Verify:_ coverage report.
- **NFR‑201 [Must]:** `README.md`, examples, and CHANGELOG updated for every feature PR.  
  _Verify:_ CI doc check; PR template gate.
- **NFR‑202 [Should]:** ADRs for architectural choices affecting APIs or performance.  
  _Verify:_ ADR directory updated.

## 6. Out of Scope (v1)
- Time‑domain solvers (FDTD/FIT), hybrid MoM/FEM, thermal/EM co‑sim  
- Massive phased‑array acceleration (beyond domain decomposition)  
- Full PCB stackup import (ODB++/Gerber) — tracked for v1.x  
- Advanced conductor roughness models — tracked for v1.x

## 7. Validation & Benchmarks
- Provide a **golden benchmark pack (≥ 20 cases)** covering: TLs (coax/microstrip), WR‑90, cavity eigenmodes/Q, patch antenna, multi‑pole filter.  
- Each benchmark records geometry hash, mesh stats, solver settings, and expected tolerances (see NFR‑100..104).  
- Nightly CI job executes a smoke subset; full pack runs on a schedule or before releases.

## 8. Traceability
- Each requirement maps to tests and sprints: see `docs/traceability.csv` with columns: `ReqID, Sprint, TestID, VerificationMethod`.  
- Example entries:
  - `FR‑010, Sprint 1, test_maxwell_mms, automated`
  - `FR‑020, Sprint 2, test_wr90_port, automated`
  - `NFR‑181, Sprint 0, scalar_smoke_macos, automated`

## 9. Acceptance & Release Criteria
- All **Must** requirements verified.  
- No open **High** severity defects in core solver/ports/PML.  
- Golden benchmark tolerances met; public examples run end‑to‑end.  
- Docs complete (README, examples, CHANGELOG); version tagged and artifacts published.
