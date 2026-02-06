# System Design Specification

**Project:** EdgeFEM — 3D FEM Electromagnetics Simulator
**Primary Modes:** Driven-modal / Driven-terminal / Eigenmode
**Outputs:** S-parameters, port impedances, fields, near/far radiation patterns, gains, efficiencies

## 1) Problem Statement & Scope

Build a production-quality 3D electromagnetic simulator based on the finite element method (FEM) for frequency-domain Maxwell’s equations in the RF–mmWave regime (≈100 kHz–110 GHz initially). Must support waveguide/coax/planar ports, lossy materials, lumped elements, PML radiation boundaries, and antenna pattern extraction. Target parity with “daily driver” HFSS workflows: filters, connectors, packages, antennas, and small arrays.

### Goals

* Accurate S-parameters (|S|, ∠S, Z/Y/ABCD) over user-defined frequency sweeps
* Robust far-field patterns (gain, directivity, efficiency, cross-pol, AR, EIRP)
* Automated, reliable meshing + hp-adaptive refinement until error/goal met
* Scales from laptop to HPC cluster; deterministic reproducibility
* Scriptable & automatable (Python API + headless CLI)

### Non-Goals (v1)

* Transient FDTD/FIT; hybrid MoM/FEM; thermal/EM co-sim; multiphysics
* Massive phased arrays with mutual coupling reduction beyond domain decomposition (add later)

## 2) Representative Use Cases

* 10-pole microstrip filter: 1–5 GHz, 3D stackup, conductor roughness
* Waveguide iris filter: X-band, 3D curved geometry, high-Q eigenmodes
* SMA/launch optimization: de-embedding to coaxial reference plane
* Antenna: patch / PIFA / helical / small array; efficiency & pattern masks
* SiP package: vias + bondwires; S-params to 40 GHz; loss & dispersion

## 3) High-Level Architecture

```
                +---------------------------+
GUI (Qt)  <-->  |  Application Layer        |  <-- Python SDK (pybind11)
CLI (headless)   +---------------------------+
                  | Model/Project Services  |  (projects, params, sweeps)
                  | Geometry & CAD Kernel   |  (import/boolean/curves)
                  | Mesh System (3D/2D)     |  (tetra/prisms, boundary mesh)
                  | EM Problem Builder      |  (ports, BCs, materials)
                  | Solver Orchestrator     |  (job graph, checkpoints)
                  | HPC Runtime             |  (MPI, GPU)
                  | Post-Processing         |  (S, fields, far-field)
                  | Results DB              |  (parquet/arrow + HDF5)
                  +---------------------------+
```

## 4) Core Numerical Method

* **Formulation:** Frequency-domain curl–curl FEM for vector E using **Nédélec edge elements** (1st–3rd order). Mixed formulation with penalty for divergence control as needed.
* **Weak form:**
  Find **E** ∈ V such that ∫(μ⁻¹ (∇×E)·(∇×v) − ω²ε E·v) dV = ∫J·v dV + boundary/port terms
* **Elements:** Curved tetrahedra (isoparametric), prisms optional for layered media; p-adaptivity (orders 1–3 v1).
* **Conductors:** Surface impedance boundary condition (SIBC) with **skin-depth** & **Huray/Groiss** roughness models (toggle).
* **Radiation:** **PML** (stretched coordinates) with automatic thickness & grading; fallback to first-order ABC for coarse passes.
* **Ports:** Wave ports (rect/wg/coax/CPW/microstrip), **driven-modal** and **driven-terminal** definitions; **lumped** ports/elements.
* **Eigenmodes:** Shift-invert for cavity/waveguide modes, Q estimation via losses.

## 5) Frequency Sweep Strategies

* **Discrete direct** (per-frequency solve) for robustness
* **Model Order Reduction (MOR):** Vector fitting / Arnoldi (Ritz vectors reuse) for smooth S(f)
* **Asymptotic/CLP** (cascaded linear predictors) for narrowband high-Q
  User chooses automatic policy: *robust*, *balanced*, or *fast* based on geometry metrics and target Q.

## 6) Linear Algebra & Solvers

* **Matrix type:** Complex, sparse, Hermitian-indefinite
* **Direct:** Multifrontal LDLᵀ/QR (e.g., MUMPS/PaStiX-like; our vendor-agnostic interface)
* **Iterative:** GMRES/FGMRES + AMG(μ⁻¹curlcurl + ω²ε) preconditioner; auxiliary-space Maxwell preconditioning
* **Parallelism:** Hybrid MPI + threads; domain decomposition (BDD/FETI-like) with coarse space
* **GPU:** Optional SpMV, triangular solves, AMG smoothers (CUDA/HIP backends; v1 focus NVIDIA)

## 7) Adaptive Refinement

* **Error estimators:** Residual-based and **goal-oriented** (S-parameter or far-field functional)
* **h-adaptivity:** Local tetra refinement with conformity (1-irregular rule)
* **p-adaptivity:** Raise polynomial order in high-curl regions; hp policy based on smoothness indicators
* **Stop criteria:** |ΔS| < ε\_S across last two refinements and/or estimator < ε\_abs; wall-clock cap

## 8) Boundary Conditions & Sources

* Perfect E/H (PEC/PMC), symmetry planes
* SIBC with σ(f), roughness, temperature coeffs
* PML regions (auto-placement from open-set faces)
* Ports: wave (eigenmode), terminal (circuit-consistent), lumped R/L/C
* Field excitations: incident plane wave (for RCS, later release)
* **De-embedding:** Port calibration to user ref planes

## 9) Materials & Libraries

* **Dielectrics:** εr(f), tanδ(f); Debye/Lorentz multi-pole models
* **Metals:** σ(f,T), roughness; option for finite thickness plating
* **Magnetics:** μr(f) for RF components (rare; supported)
* **Stackups:** PCB layer manager (copper/FR-4/Rogers), substrate anisotropy
* **Library:** Versioned material DB (JSON + checksum), import from CSV

## 10) Geometry & Meshing

* **CAD import:** STEP, IGES, Parasolid; PCB: ODB++, Gerber+stackup
* **Ops:** Booleans, fillets, defeaturing (slivers, tiny edges), parameterization
* **Meshing:**

  * **Volume:** Delaunay/front-advancing tetra; curvature & size fields; boundary layers for conductors
  * **Surface:** Curvature-adapted triangles; port cross-section mesher (eigen mesh quality > target)
  * Quality metrics: min dihedral, aspect ratio, Jacobian > 0.3; smart sliver removal

## 11) Post-Processing & Observables

* **S-parameters:** Magnitude/phase, Smith, group delay, Q, pass/stop masks
* **Fields:** |E|, |H|, J, power flow; cut-planes, streamlines; field probes
* **Loss accounting:** dielectric/conductor/PML; radiation & total efficiency
* **Antenna patterns:**

  * Near-to-far via equivalent currents on a Huygens surface or Stratton–Chu
  * 2D/3D patterns; RHCP/LHCP, AR, cross-pol; EIRP; gain/directivity; beamwidth; sidelobe levels
  * **Pattern integration** against regulatory masks (basic v1)
* **Reports:** Templates (filter, interconnect, antenna); multi-design comparators
* **Data export:** Touchstone (.sNp), CSV, HDF5/Arrow for fields & meshes, .pat JSON

## 12) Results & Data Model

* **Project store:**

  * Geometry: CAD + param set (JSON)
  * Mesh: HDF5 (nodes, edges, elems, p-orders, regions)
  * Solver state: checkpoints per frequency (for MOR reuse)
  * Results DB: Parquet/Arrow + catalog metadata for fast slicing (freq, port pair, probe)
* **Reproducibility:** Hash of CAD+mesh+settings; run manifest; semantic versioning

## 13) HPC Runtime & Job Management

* **Scheduler integrations:** SLURM, PBS, local pool
* **Job graph:** Preprocess → Mesh → (freq nodes) Solve → Post → Reduce/MOR → Reports
* **Checkpointing:** After factorization / Krylov basis; fault-tolerant restart
* **Licensing server (optional):** Token pool per concurrent solve

## 14) APIs, Scripting, and Automation

* **Python SDK:**

  * `project = em.Project()`
  * `g = project.geometry.import_step("cavity.step")`
  * `project.materials.assign(g.faces, "Copper", roughness="huray")`
  * `port = project.ports.wave(face=port_face, modes=1)`
  * `sweep = project.sweep.linear(1e9, 5e9, 401, policy="balanced")`
  * `res = project.solve(hp_adapt=True)`
  * `s = res.sparams(); pat = res.farfield(f=3e9)`
* **CLI:** `edgefem solve project.veproj --hpadapt --scheduler slurm`
* **Plugins:** C++/Python hooks for custom materials, post-processors

## 15) Verification & Validation

* **Unit tests:** Element matrices, BCs, port eigenmodes, PML absorption
* **Method of Manufactured Solutions (MMS)** for regression
* **Benchmarks:** Rect/coax lines (characteristic Z), WG steps, cavity eigenfrequencies/Q, patch antenna gain vs. measurement
* **Tolerances:**

  * S-param |ΔS| < 0.01 for mesh refinement step-to-step
  * Eigenfreq < 0.5% error; Q within 5–10% depending on loss model
  * Pattern total radiated power balance within 1%

## 16) Performance Targets (v1 hardware assumptions)

* 8-core workstation: 2–5M dof, single frequency in <30 min (direct), <10 min (iterative good precond)
* 2×GPU node: speedup 2–4× vs CPU for iterative solves
* Cluster 64 ranks: strong scaling efficiency ≥70% to 50M dof
* Memory headroom: 16 bytes/nnz target; peak < 2.5× matrix for direct solves

## 17) UX & Workflow

* **GUI:** HFSS-style tree (Geometry/Setup/Excitations/Analysis/Results), property panels, interactive plots (Bode, Smith, polar, 3D lobes)
* **Wizards:** Transmission line, waveguide, antenna templates
* **Design variables & param sweeps:** DOE/Opt (Sobol, Nelder–Mead; hooks for external optimizers)
* **Diagnostics:** Mesh-quality heatmaps, PML grading plots, port-mode field preview

## 18) Configuration & Deployment

* **Languages:** C++20 (core), CUDA/HIP (accel), Python 3.11 (SDK/GUI glue), Qt 6
* **Build:** CMake; conda env for Python bits
* **Platforms:** Linux (primary), Windows (secondary), macOS (dev only, no GPU v1)
* **Packaging:** Installer + runtime redistributables; container images for HPC

## 19) Observability

* **Logging:** Structured (spdlog) with run IDs; levels per subsystem
* **Telemetry (opt-in):** Anon performance counters, crash dumps
* **Profiling:** Timers around mesh/assemble/solve; flamegraphs in GUI

## 20) Security & IP

* Encrypted project archives (optional)
* Checksummed material libs; signed plugins
* No external net access during solves (HPC safe)

## 21) Roadmap (abridged)

* **v1.0:** Driven-modal/terminal, ports, PML, S-params, far-field, hp-adapt, direct + GMRES+AMG, Python API, SLURM
* **v1.1:** Better MOR, GPU AMG, conductor roughness models, batch optimization
* **v1.2:** Eigenmode Q improvements, PCB stackup import, pattern mask compliance
* **v1.3+:** Hybrid FEM-MoM, RCS, thermal coupling, array domain decomposition

## 22) Risks & Mitigations

* **PML instability at grazing incidence:** auto-tuning & validation suite
* **AMG for curl-curl is tricky:** auxiliary-space + robust defaults; fall back to direct
* **Mesh curvature on tight radii:** curved elements + geometry tolerance cleanup
* **User error in port setup:** field previews, port-QC checks (orthogonality/leakage)

## 23) Acceptance Criteria (v1 release)

* Golden set of 20 benchmarks passing tolerances (S, eigen, patterns)
* End-to-end regression time < 2 hrs CI on standard runner
* Docs: user guide, Python examples, validation notebook pack

---

### Appendix A — Key Data Structures

* `Mesh`: CSR adjacency + element arrays (conn, region id, p-order, curved map)
* `ElementOps`: shape funcs, curl/curl stencils, numerical quadrature
* `BC/Port`: typed struct with references to boundary facets & eigen basis
* `SystemMatrix`: block operators; lazy assembly; views for preconditioners
* `ResultSet`: Arrow/Parquet tables keyed by `run_id, freq, port_i, port_j`

### Appendix B — Default Parameters

* PML thickness = max(λ/4, 5× elem size), grading = polynomial(3), target refl ≤ −40 dB
* h-refine if local estimator > 0.7×max; p-refine if smoothness > threshold
* AMG: V-cycle depth 3, Chebyshev smoother; restart GMRES(60)

---

If you want, I can turn this into:

* a **requirements doc (PRD)** for stakeholders,
* a **tasks + milestones plan** for engineering,
* or **Python stubs** for the SDK to kickstart prototyping.

Any special constraints (target bands, must-have CAD formats, GPU availability, license model)?
