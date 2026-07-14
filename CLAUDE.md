# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EdgeFEM is a 3D finite-element electromagnetics (FEM) simulator for **waveguide, patch antenna, and periodic unit cell RF modeling** (validated ~5–26 GHz). It computes S-parameters, fields, and radiation patterns using lowest-order Nédélec (edge) elements, with single-axis Bloch periodic BCs and a graded absorbing layer for open boundaries. Execution is single-threaded, single-node; direct-solver memory limits apply.

Development tooling includes small Python CLIs under `tools/` (mesh QA, simulation scaffolding, pattern plotting); see AGENTS.md. There is no MCP server.

## Ecosystem Role

EdgeFEM serves as the **full-wave FEM engine** in a multi-package RF modeling ecosystem:

```
┌─────────────────────────────────────────────────────────────────┐
│                     Application Layer                           │
│  ┌─────────────────────┐    ┌─────────────────────────────────┐ │
│  │ Phased-Array-       │    │ RF Metasurface/Metamaterial     │ │
│  │ Antenna-Model       │    │ Package (planned)               │ │
│  │ - Array factor      │    │ - Unit cell optimization        │ │
│  │ - Beamforming       │    │ - Homogenization                │ │
│  │ - Scan analysis     │    │ - Surface impedance extraction  │ │
│  └─────────┬───────────┘    └───────────────┬─────────────────┘ │
│            │                                │                   │
│            ▼                                ▼                   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                      EdgeFEM                             │   │
│  │  Full-wave 3D FEM engine providing:                       │   │
│  │  - Unit cell S-parameters (periodic BCs)                  │   │
│  │  - Embedded element patterns                              │   │
│  │  - Mutual coupling coefficients                           │   │
│  │  - Near-field distributions                               │   │
│  └──────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

### Interoperability Interfaces

| Interface | Format | Use Case |
|-----------|--------|----------|
| Python SDK (`pyedgefem`) | NumPy arrays | Direct integration with array/metasurface packages |
| Touchstone (.sNp) | IEEE standard | S-parameter exchange with any RF tool |
| JSON | Structured data | Element patterns, coupling matrices, geometry |
| Gmsh v2 | Mesh format | Geometry import from CAD pipelines |

### Integration with Phased-Array-Antenna-Model

EdgeFEM provides element-level data that feeds into [Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model):
- **Element patterns**: Far-field radiation from unit cell simulation → array factor weighting
- **Mutual coupling**: Full-wave S-matrix between adjacent elements → impairment modeling
- **Scan impedance**: Port impedance vs. scan angle → active VSWR prediction

### Integration with Metasurface Package (planned)

EdgeFEM will provide unit cell characterization for metasurface/metamaterial design:
- **Reflection/transmission coefficients**: Floquet-port S-parameters under oblique incidence
- **Surface impedance**: Extracted Zs for equivalent surface models
- **Dispersion curves**: Frequency-dependent effective parameters (ε_eff, μ_eff)

## Build Commands

```bash
# macOS prerequisites
brew install cmake ninja eigen

# Release build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# With Python SDK
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DEDGEFEM_PYTHON=ON
cmake --build build -j
```

## Test Commands

```bash
# Smoke tests (fast)
ctest --test-dir build -L smoke -j

# All tests
ctest --test-dir build -j

# Specific test labels
ctest --test-dir build -L pml    # PML absorption
ctest --test-dir build -L abc    # ABC boundary

# Single test binary
./build/tests/test_maxwell
```

## Lint Commands

```bash
# Format and static analysis
bash tools/lint.sh build

# clang-tidy only
clang-tidy -p build -warnings-as-errors='*' $(git ls-files 'src/*.cpp')
```

## Demo Commands

```bash
# Scalar Helmholtz cavity
./build/src/edgefem_scalar_demo examples/cube_cavity.msh 1

# Waveguide S-parameters
./build/src/edgefem_waveguide_demo
python3 examples/plot_waveguide_sparams.py
```

## Architecture

### Core Modules (`include/edgefem/`)

| Module | Purpose |
|--------|---------|
| `mesh.hpp` | Tet4/Tri3 mesh with global edge indexing, loaded from Gmsh v2 |
| `maxwell.hpp` | Maxwell curl-curl FEM assembly with PML and ABC |
| `edge_basis.hpp` | Nédélec (Whitney) edge element basis functions |
| `ports/wave_port.hpp` | Wave port definitions with 2D eigenmode solver |
| `solver.hpp` | BiCGSTAB+ILUT solver with auto-fallback to SparseLU |
| `sweep.hpp` | Frequency sweep with K/M separation and factorization reuse |
| `ports/lumped_port.hpp` | Lumped port with surface-integral weight computation |
| `io/touchstone.hpp` | S-parameter Touchstone export |
| `post/ntf.hpp` | Near-to-far field via Huygens surface |
| `mor/vector_fit.hpp` | Rational model fitting for MOR |

### Data Flow

1. **Mesh** → `load_gmsh_v2()` parses Gmsh file, builds edge connectivity
2. **BC** → `build_edge_pec()` tags PEC edges from physical groups
3. **Ports** → `build_wave_port_from_eigenvector()` for wave ports, `build_lumped_port()` for antennas
4. **Assembly** → `assemble_maxwell()` or `assemble_maxwell_km()` (K/M separation for sweeps)
5. **Solve** → `solve_linear()` via BiCGSTAB+ILUT with auto-fallback to SparseLU
6. **Post** → `calculate_sparams_eigenmode()` for wave ports, `calculate_sparams()` for lumped ports

### Key Types

- `VecC = Eigen::VectorXcd` — complex solution vectors
- `SpMatC = Eigen::SparseMatrix<std::complex<double>>` — sparse system matrices
- Edge orientation stored per element: `edge_orient[i]` ∈ {+1, −1}

## Developer Tooling

Python CLIs under `tools/` (see AGENTS.md):
- `tools/agents/mesh_qa.py` — Gmsh v2 mesh quality checks (volumes, aspect ratios), JSON output; runs in CI
- `tools/new_simulation.py` — scaffolds a `.geo` + `run_simulation.py` from templates
- `tools/plot_pattern.py` — polar pattern plots from CSV

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `EDGEFEM_BUILD_SCALAR` | ON | Scalar Helmholtz path |
| `EDGEFEM_BUILD_VECTOR` | ON | Maxwell vector path |
| `EDGEFEM_PYTHON` | OFF | pybind11 Python bindings |

## Code Conventions

- All code in `namespace edgefem`
- Public headers in `include/edgefem/`, implementations in `src/`
- Tests use simple `assert()` macros (no external framework)
- Edge canonical form: `(n0 < n1)` stored in `mesh.edges`

## Known Issues & Critical Fixes

### Wave Port: Assembled Surface Mass Matrix + 2D Discrete Port Mode (July 2026)

**Supersedes the 0.5 ABC scale factor.** `calculate_sparams_eigenmode()` now:
1. Applies the assembled port surface mass matrix: `A += jβ·M_s` per port
   (`assemble_port_surface_mass()` in `src/ports/wave_port.cpp`), replacing
   the diagonal-lumped ABC that required the empirical `port_abc_scale = 0.5`.
2. Uses the 2D discrete port-face eigenmode as port weights
   (`build_wave_port_2d()` / `solve_port_mode_2d()`): a dense generalized
   EVP (K_s, M_s) on the port edges. `mode.kc` is replaced by the discrete
   cutoff so β is consistent with the profile.
3. Normalizes and extracts in the M_s inner product (`e^H M_s e = 1`,
   `V_j = e_j^H M_s x`).

Measured on WR-90 (7-12 GHz): |S11| 0.003-0.053 (was 0.04-0.14),
|S21| ≥ 0.997, ABC-scale sweep optimum exactly 1.0. `port_abc_scale` now
defaults to 1.0 and is a diagnostic knob only — do not tune it.

**Guidelines**: build ports with `build_wave_port_2d(mesh, tag, mode,
pec_edges, kc_sq_target)`. The older `compute_te_eigenvector` +
`build_wave_port_from_eigenvector` path (3D cavity eigensolve) is slower
and its profile carries ~10% non-modal contamination; avoid for new code.

**Critical convention**: the mesh loader (`build_edges()` in
`src/mesh_gmsh.cpp`) fills Tri3 `Element::edges` in local order
**(0,1), (1,2), (2,0)**. Any surface assembly must use that table. A
mismatched table `{{0,1},{0,2},{1,2}}` in `lumped_port.cpp` silently
scattered surface integrals to wrong edges until July 2026 — validate new
surface code with the constant-field zero-curl test in
`tests/test_triangle_mass_matrix.cpp`.

### Wave Port Formulation: 2D vs 3D Mode Coupling (Fixed Feb 2026, superseded July 2026)

**Problem**: S-parameter calculations produced non-physical results (|S11| > 1, |S21| ≈ 0) for simple waveguide structures that should have |S11| ≈ 0, |S21| ≈ 1.

**Root Cause**: 2D analytical port weights had only **3% correlation** with the 3D FEM eigenvector at port edges.

| Component | Expected | Actual | Issue |
|-----------|----------|--------|-------|
| Analytical weights Y/X edge ratio | ~12.7 | 0.95 | Analytical formula distributes weight evenly |
| 3D FEM eigenvector Y/X edge ratio | 12.7 | 12.7 | Correct (y-edges dominate for TE10) |
| Weight correlation |<w_ana, v_3D>| | ~1.0 | 0.03 | Fundamental mismatch |

**Symptoms indicating this issue**:
- |S11| > 1 (passivity violation)
- |S21| ≈ 0 when propagation expected
- Solution field decays along waveguide instead of propagating
- Port loading (ww^H/Z0) << 1% of system matrix norm

**Implemented Fixes**:

1. **Eigenvector-based port weights** (`src/ports/wave_port.cpp`):
   ```cpp
   // Use 3D FEM eigenvector instead of analytical formula
   WavePort build_wave_port_from_eigenvector(mesh, surface, eigenvector, mode, pec_edges);
   Eigen::VectorXd compute_te_eigenvector(mesh, pec_edges, target_kc_sq);
   ```

2. **Direct eigenmode S-parameter calculation** (`src/assemble_maxwell.cpp`):
   ```cpp
   // Use eigenmode-based S-parameter extraction with optimal ABC scaling
   MaxwellParams p;
   p.omega = 2 * M_PI * freq;
   p.port_abc_scale = 0.5;  // Optimal scaling factor
   Eigen::MatrixXcd S = calculate_sparams_eigenmode(mesh, p, bc, ports);
   ```

   The key insight: the optimal ABC coefficient is 0.5×jβ (not 1.0×jβ) due to the
   diagonal approximation of the surface mass matrix in the FEM discretization.

**Guidelines for future port development**:
- Use `calculate_sparams_eigenmode()` for accurate S-parameters (>99% transmission accuracy)
- Set `MaxwellParams::port_abc_scale = 0.5` (optimal value from tuning study)
- Always verify port weight correlation with 3D FEM eigenvector: |<w, v>| should be > 0.9
- Check edge direction bias: for TE10, y-directed edges should dominate
- Test passivity: |S11|² + |S21|² ≤ 1.01
- Verify S21 phase matches analytical: phase(S21) = -β×L (within ~10° for typical meshes)

**Related test files**:
- `tests/test_eigenvector_waveport.cpp` — validates eigenvector-based approach
- `tests/debug_mode_alignment.cpp` — diagnoses weight/eigenvector correlation
- `tests/debug_3d_vs_2d_pattern.cpp` — compares edge direction distributions

### Mode Normalization for Unit Power Flow (Fixed Feb 2026)

**Problem**: TE10 mode amplitude was 87.5× too small due to missing factor in normalization formula.

**Root Cause**: Original formula `A² = 4·kc⁴/(ω·μ·β·π²·b)` was missing factor `a` (waveguide width).

**Fix**: Corrected to `A² = 4·kc⁴·a/(ω·μ·β·π²·b)` in `src/ports/wave_port.cpp`.

**Derivation**: For TE10 mode with Hz = A·cos(πx/a), unit power flow requires:
```
P = (1/2)·Re∫∫ E×H* · dS = 1
→ A² = 4·kc⁴·a / (ω·μ·β·π²·b)
```

### Lumped Port Passivity & Weight Normalization (Fixed Mar 2026)

**Problem**: Stacked patch antenna showed |S11| = 5.72 (severe passivity violation). The `simulate()` method used a hand-rolled extraction with magnitude normalization and passivity clamping that masked the real issue.

**Root Cause**: `build_lumped_port()` used a crude `edge_vec.dot(e_dir)` projection for port weights, and `normalize_port_weights()` was never called before `calculate_sparams()`.

**Fixes**:
1. **Surface-integral weight computation** (`src/ports/lumped_port.cpp`): Whitney basis mass matrix integration `w_i = ∫(N_i · e_dir) dS` replaces raw edge projection. Set `LumpedPortConfig::weight_mode = SurfaceIntegral` (now the default).
2. **Standard S-parameter path** (`python/edgefem/designs/stacked_patch.py`): `simulate()` now calls `normalize_port_weights()` → `calculate_sparams()` instead of hand-rolled extraction.
3. **ILUT preconditioner** (`src/solver.cpp`): BiCGSTAB now uses `Eigen::IncompleteLUT` by default. Auto-fallback to SparseLU on convergence failure.
4. **Frequency sweep efficiency** (`src/sweep.cpp`): `frequency_sweep()` pre-assembles K and M separately, reuses symbolic factorization across frequencies.

**Guidelines for lumped port development**:
- Always call `normalize_port_weights()` before `calculate_sparams()` for lumped ports
- Use `SurfaceIntegral` weight mode (default) for mesh-convergent weights
- Use `convergence_check()` on design classes to verify mesh adequacy
- Passivity violations (|S11| > 1) on coarse meshes are expected; refine mesh rather than clamping

### Waveguide Port/Eigenmode Path (Fixed Mar 2026)

**Problem**: `RectWaveguideDesign` produced |S21| = 0.014 instead of |S21| > 0.95.

**Root Cause**: `_setup_ports()` used `build_wave_port()` (analytical weights) but passed to `calculate_sparams_eigenmode()` which expects `build_wave_port_from_eigenvector()` (FEM eigenvector weights with `||w||² = sqrt(Z0)` normalization).

**Fix**: `_setup_ports()` now uses `compute_te_eigenvector()` → `build_wave_port_from_eigenvector()`, matching the proven C++ test path.

## Reference Documents

- `system_design_spec.md` — numerical methods, architecture decisions
- `project_plan.md` — sprint breakdown, acceptance criteria
- `AGENTS.md` — automation agent design, MCP wiring
