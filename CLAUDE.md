# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VectorEM is a 3D finite-element electromagnetics (FEM) simulator specialized for **metasurface, metamaterial, and phased array unit cell RF modeling** (100 kHz–110 GHz). It computes S-parameters, fields, and radiation patterns using Nédélec (edge) elements with periodic boundary conditions, PML for open boundaries, and scales from laptops to HPC clusters.

The project is MCP (Model Context Protocol) compatible, enabling LLM-augmented automation agents for mesh QA, numerical validation, and documentation sync.

## Ecosystem Role

VectorEM serves as the **full-wave FEM engine** in a multi-package RF modeling ecosystem:

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
│  │                      VectorEM                             │   │
│  │  Full-wave 3D FEM engine providing:                       │   │
│  │  - Unit cell S-parameters (Floquet ports)                 │   │
│  │  - Embedded element patterns                              │   │
│  │  - Mutual coupling coefficients                           │   │
│  │  - Near-field distributions                               │   │
│  └──────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

### Interoperability Interfaces

| Interface | Format | Use Case |
|-----------|--------|----------|
| Python SDK (`pyvectorem`) | NumPy arrays | Direct integration with array/metasurface packages |
| Touchstone (.sNp) | IEEE standard | S-parameter exchange with any RF tool |
| JSON | Structured data | Element patterns, coupling matrices, geometry |
| Gmsh v2 | Mesh format | Geometry import from CAD pipelines |

### Integration with Phased-Array-Antenna-Model

VectorEM provides element-level data that feeds into [Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model):
- **Element patterns**: Far-field radiation from unit cell simulation → array factor weighting
- **Mutual coupling**: Full-wave S-matrix between adjacent elements → impairment modeling
- **Scan impedance**: Port impedance vs. scan angle → active VSWR prediction

### Integration with Metasurface Package (planned)

VectorEM will provide unit cell characterization for metasurface/metamaterial design:
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
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DVECTOREM_PYTHON=ON
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
./build/src/vectorem_scalar_demo examples/cube_cavity.msh 1

# Waveguide S-parameters
./build/src/vectorem_waveguide_demo
python3 examples/plot_waveguide_sparams.py
```

## Architecture

### Core Modules (`include/vectorem/`)

| Module | Purpose |
|--------|---------|
| `mesh.hpp` | Tet4/Tri3 mesh with global edge indexing, loaded from Gmsh v2 |
| `maxwell.hpp` | Maxwell curl-curl FEM assembly with PML and ABC |
| `edge_basis.hpp` | Nédélec (Whitney) edge element basis functions |
| `ports/wave_port.hpp` | Wave port definitions with 2D eigenmode solver |
| `solver.hpp` | BiCGSTAB + ILUT linear solver interface |
| `sweep.hpp` | Frequency sweep orchestration |
| `io/touchstone.hpp` | S-parameter Touchstone export |
| `post/ntf.hpp` | Near-to-far field via Huygens surface |
| `mor/vector_fit.hpp` | Rational model fitting for MOR |

### Data Flow

1. **Mesh** → `load_gmsh_v2()` parses Gmsh file, builds edge connectivity
2. **BC** → `build_edge_pec()` tags PEC edges from physical groups
3. **Ports** → `build_wave_port()` extracts surface mesh, solves 2D eigenmode
4. **Assembly** → `assemble_maxwell()` builds curl-curl + mass matrices with PML/ABC
5. **Solve** → `solve_linear()` via BiCGSTAB with ILUT preconditioner
6. **Post** → `calculate_sparams()` extracts S-matrix, exports to Touchstone

### Key Types

- `VecC = Eigen::VectorXcd` — complex solution vectors
- `SpMatC = Eigen::SparseMatrix<std::complex<double>>` — sparse system matrices
- Edge orientation stored per element: `edge_orient[i]` ∈ {+1, −1}

## MCP Integration

Automation agents in `tools/agents/` follow the MCP server pattern:
- `mesh_qa.py` — validates mesh quality, outputs JSON diagnostics
- `release_notes.py` — generates changelog from merged PRs
- Designed for LLM-augmented workflows with structured output

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `VECTOREM_BUILD_SCALAR` | ON | Scalar Helmholtz path |
| `VECTOREM_BUILD_VECTOR` | ON | Maxwell vector path |
| `VECTOREM_PYTHON` | OFF | pybind11 Python bindings |

## Code Conventions

- All code in `namespace vectorem`
- Public headers in `include/vectorem/`, implementations in `src/`
- Tests use simple `assert()` macros (no external framework)
- Edge canonical form: `(n0 < n1)` stored in `mesh.edges`

## Known Issues & Critical Fixes

### Wave Port Formulation: 2D vs 3D Mode Coupling (Fixed Feb 2026)

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

## Reference Documents

- `system_design_spec.md` — numerical methods, architecture decisions
- `project_plan.md` — sprint breakdown, acceptance criteria
- `AGENTS.md` — automation agent design, MCP wiring
