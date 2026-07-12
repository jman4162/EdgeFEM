# EdgeFEM

**3D Finite-Element Electromagnetics Simulator for RF and mmWave**

[![PyPI version](https://img.shields.io/pypi/v/edgefem.svg)](https://pypi.org/project/edgefem/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![CI](https://github.com/jman4162/EdgeFEM/actions/workflows/ci.yml/badge.svg)](https://github.com/jman4162/EdgeFEM/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://jman4162.github.io/EdgeFEM/)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

EdgeFEM is an open-source full-wave frequency-domain FEM solver for RF and microwave structures: hollow waveguides, probe-fed patch antennas, and periodic unit cells. It computes S-parameters, electromagnetic fields, and radiation patterns using lowest-order Nédélec (edge) tetrahedral elements. Accuracy is validated against analytical benchmarks over roughly 5–26 GHz (see [Validation](docs/validation.md)).

## Features

| Feature | Description |
|---------|-------------|
| **S-parameters** | Wave ports with 3D-eigenvector-based extraction; lumped ports for antenna feeds |
| **Radiation Patterns** | 3D far-field via Stratton-Chu integration over a Huygens surface |
| **Periodic Structures** | Phase-shifted (Bloch) periodic BCs along one lattice axis, for small unit cells |
| **Open Boundaries** | First-order ABC and a graded absorbing layer (see Limitations) |
| **Dispersive Materials** | Debye, Lorentz, Drude, and DrudeLorentz models |
| **Python SDK** | Full scriptable API with NumPy integration |
| **Design Classes** | `RectWaveguideDesign`, `PatchAntennaDesign`, `UnitCellDesign` |
| **Plotting Module** | 14 visualization functions |
| **Standard Exports** | Touchstone (.sNp), VTK, CSV |

## Quick Start

### Installation

```bash
# macOS
brew install cmake ninja eigen gmsh

# Clone and build
git clone https://github.com/jman4162/EdgeFEM.git
cd EdgeFEM
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DEDGEFEM_PYTHON=ON
cmake --build build -j

# Set Python path
export PYTHONPATH="$PWD/build/python:$PWD/python:$PYTHONPATH"
```

### 10-Line Example: Waveguide S-Parameters

```python
from edgefem.designs import RectWaveguideDesign
import numpy as np

# Create WR-90 waveguide (X-band)
wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
wg.generate_mesh(density=10)

# Compute S-parameters at 10 GHz
S = wg.sparams_at_freq(10e9)
print(f"|S11| = {abs(S[0,0]):.4f}")  # ~0.09 (see validation report)
print(f"|S21| = {abs(S[1,0]):.4f}")  # ~0.99 (matched transmission)
```

### Patch Antenna Example (Analytical Estimates)

```python
from edgefem.designs import PatchAntennaDesign

# 2.4 GHz WiFi patch on FR-4
# Note: input_impedance() and radiation_pattern() use analytical
# approximations (transmission-line / cavity model). For full-wave FEM
# results use patch.simulate(freq) after generate_mesh().
patch = PatchAntennaDesign(
    patch_length=29e-3, patch_width=38e-3,
    substrate_height=1.6e-3, substrate_eps_r=4.4
)
patch.set_probe_feed(y_offset=-5e-3)
patch.generate_mesh(density=15)

# Analyze (analytical model — suitable for initial design)
Z_in = patch.input_impedance(2.4e9)
pattern = patch.radiation_pattern(2.4e9)
D = patch.directivity(2.4e9)
print(f"Directivity: {10*np.log10(D):.1f} dBi")
```

### Unit Cell / Metasurface Example

```python
from edgefem.designs import UnitCellDesign

# Periodic unit cell
cell = UnitCellDesign(period_x=5e-3, period_y=5e-3,
                       substrate_height=0.5e-3, substrate_eps_r=3.5)
cell.add_patch(width=4e-3, length=4e-3)
cell.generate_mesh(density=20)

# Reflection coefficient at normal incidence
R, T = cell.reflection_transmission(10e9, theta=0, phi=0)
Zs = cell.surface_impedance(10e9)
```

Note: the periodic solve currently applies the Bloch phase along the x lattice
axis only (y walls are natural boundaries), so results are most reliable at
normal incidence. The port is the cell's dominant waveguide mode, not a
Floquet harmonic. See Current Limitations.

## Documentation

- **[Installation Guide](docs/installation.md)** - Build from source
- **[Quick Start](docs/quickstart.md)** - First simulation in 5 minutes
- **[API Reference](docs/api/)** - Complete function documentation
- **[Tutorials](tutorials/)** - Jupyter notebook examples
- **[Validation Report](docs/validation.md)** - Accuracy benchmarks

### Jupyter Notebook Tutorials

| Tutorial | Description |
|----------|-------------|
| [01_rectangular_waveguide.ipynb](tutorials/01_rectangular_waveguide.ipynb) | WR-90 S-parameters |
| [02_patch_antenna.ipynb](tutorials/02_patch_antenna.ipynb) | 2.4 GHz patch design |
| [03_unit_cell.ipynb](tutorials/03_unit_cell.ipynb) | Metasurface characterization |
| [04_phased_array.ipynb](tutorials/04_phased_array.ipynb) | Active impedance vs scan |

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                      Python SDK (pyedgefem)                     │
│  ┌─────────────────────┐  ┌────────────────┐  ┌──────────────┐  │
│  │   Design Classes    │  │  Plots Module  │  │  Core API    │  │
│  │ RectWaveguideDesign │  │ plot_pattern   │  │ load_gmsh    │  │
│  │ PatchAntennaDesign  │  │ plot_sparams   │  │ assemble     │  │
│  │ UnitCellDesign      │  │ plot_coupling  │  │ solve        │  │
│  └─────────────────────┘  └────────────────┘  └──────────────┘  │
└─────────────────────────────┬───────────────────────────────────┘
                              │ pybind11
┌─────────────────────────────┴───────────────────────────────────┐
│                          C++ Core                                │
│  ┌──────────┐ ┌──────────┐ ┌────────┐ ┌──────────┐ ┌─────────┐  │
│  │   Mesh   │ │  Maxwell │ │ Ports  │ │  Solver  │ │  Post   │  │
│  │  (Gmsh)  │ │ Assembly │ │ (Wave/ │ │(BiCGSTAB)│ │ (NTF,   │  │
│  │          │ │(ABC/absr)│ │ Lumped)│ │          │ │  VTK)   │  │
│  └──────────┘ └──────────┘ └────────┘ └──────────┘ └─────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

## Validation

Validated against analytical benchmarks (WR-90/WR-42 waveguide dispersion, cavity eigenmodes, lossy-dielectric attenuation). Representative WR-90 results at ~10 elements/wavelength:

| Metric | Typical (WR-90, 7-12 GHz) | Test pass criterion |
|--------|---------------------------|---------------------|
| \|S21\| | 0.987-0.998 | > 0.90 |
| \|S11\| | 0.04-0.14 | < 0.15 |
| Phase error vs. -βL | 3-13° | < 15° |
| Passivity \|S11\|² + \|S21\|² | ≈ 1.0 | ≤ 1.05 |

Numbers are from `tests/benchmark_wr90.cpp`, run in CI. See [Validation Report](docs/validation.md) for full results and limitations.

## Ecosystem Integration

EdgeFEM is the full-wave FEM engine in a multi-package RF modeling ecosystem:

```
┌───────────────────────────────────────────────────────────────┐
│                    Application Layer                          │
│  ┌─────────────────────┐  ┌─────────────────────────────────┐ │
│  │ Phased-Array-       │  │ RF Metasurface Package          │ │
│  │ Antenna-Model       │  │ (planned)                       │ │
│  └──────────┬──────────┘  └────────────────┬────────────────┘ │
└─────────────┼──────────────────────────────┼──────────────────┘
              │                              │
              ▼                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                          EdgeFEM                                │
│  Full-wave 3D FEM providing: S-params, patterns, coupling       │
└─────────────────────────────────────────────────────────────────┘
```

- **[Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model)**: Uses EdgeFEM element patterns and coupling for array synthesis

## Build Options

```bash
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DEDGEFEM_PYTHON=ON         # Python bindings
```

| Option | Default | Description |
|--------|---------|-------------|
| `EDGEFEM_BUILD_SCALAR` | ON | Scalar Helmholtz solver |
| `EDGEFEM_BUILD_VECTOR` | ON | Maxwell vector solver |
| `EDGEFEM_PYTHON` | OFF | Python bindings (pybind11) |

## Run Tests

```bash
ctest --test-dir build -j        # All tests
ctest --test-dir build -L smoke  # Quick smoke tests
ctest --test-dir build -L pml    # PML absorption tests
```

## Plotting Module

```python
import edgefem.plots as vp

# S-parameters vs frequency
vp.plot_sparams_vs_freq(freqs, S_list, params=['S11', 'S21'])

# Radiation patterns
vp.plot_pattern_3d(pattern, interactive=True)
vp.plot_pattern_3d_cuts(pattern, title="E/H-plane")

# Coupling analysis
vp.plot_coupling_matrix(S, title="4×4 Array")
vp.plot_active_impedance_scan(angles, Z_active)
```

## Export Formats

| Format | Function | Use Case |
|--------|----------|----------|
| Touchstone (.sNp) | `export_touchstone()` | RF circuit simulators |
| VTK (.vtu) | `export_fields_vtk()` | ParaView visualization |
| CSV | `export_pattern_csv()` | Custom post-processing |

## Citation

If you use EdgeFEM in your research, please cite:

```bibtex
@software{hodge2026edgefem,
  author       = {Hodge, John},
  title        = {{EdgeFEM}: A 3D Finite-Element Electromagnetics Simulator for RF and mmWave},
  year         = {2026},
  month        = feb,
  version      = {1.0.0},
  publisher    = {GitHub},
  url          = {https://github.com/jman4162/EdgeFEM},
  license      = {MIT}
}
```

## Current Limitations

EdgeFEM v1.0 targets small-to-medium waveguide, patch antenna, and unit cell problems. Know these limits before relying on results:

| Feature | Status | Notes |
|---------|--------|-------|
| True tensor (coordinate-stretched) PML | Planned v1.1 | Current absorber applies a scalar average of the stretch tensor: a graded lossy layer, not reflectionless at oblique incidence |
| Wave port formulation | First-order ABC with modal injection | Diagonal (lumped) surface operator with tuned 0.5 scale factor; ~5% S-parameter error typical; accuracy degrades near cutoff |
| Floquet ports / harmonic expansion | Planned | Periodic BC applies a Bloch phase along one lattice axis only; the constraint elimination is dense, limiting problem size |
| Inhomogeneous / TEM ports (coax, microstrip) | Not supported | 2D port eigensolver is scalar (hollow homogeneous guides only) |
| Higher-order elements (Tet10) | Planned v1.1 | Lowest-order Whitney elements; use finer mesh |
| Adaptive mesh refinement | Planned v1.1 | Manual mesh refinement in Gmsh |
| Higher-order modes (TE20, TM) | Planned v1.1 | Single-mode analysis only |
| Plane wave excitation (RCS) | Planned v1.2 | — |
| Anisotropic materials | Not planned | Isotropic (complex scalar) only |
| Parallelism (threads/MPI/GPU) | Planned | Single-threaded execution; direct-solver memory limits apply |
| Low-frequency stabilization | Not implemented | Expect conditioning problems well below ~1 GHz |
| Time-domain solver | Not planned | Frequency-domain only |
| Gmsh mesh input | v2 ASCII only | Export with `gmsh -format msh2`; v4 files are not readable |

For large problems (>1M DOF), adaptive refinement, or reference-plane modal ports, use commercial tools like HFSS or CST.

## Roadmap

### v1.0 (Current)
- Dispersive materials (Debye, Lorentz, Drude, DrudeLorentz models)
- pip-installable Python package (`pip install edgefem`)
- GitHub Actions CI/CD with automated wheel builds

### v1.1 (Planned)
- True tensor (coordinate-stretched) PML
- Port surface mass matrix (replaces the tuned ABC scale factor)
- Sparse periodic constraint elimination; dual-axis periodicity
- Higher-order waveguide modes (TE20, TE01, TM modes)
- Higher-order Nédélec elements (Tet10)
- OpenMP parallelization
- Adaptive mesh refinement (h-adaptivity)

### v1.2 (Planned)
- Floquet-harmonic ports for oblique-incidence unit cell analysis
- Plane wave excitation for RCS/scattering
- Surface impedance BC (SIBC)
- Direct solver option (MUMPS)

### v2.0 (Future)
- MPI distributed solving
- GPU acceleration
- Optimization framework
- GUI application

## Contributing

- Use feature branches; PRs require green CI and updated docs
- clang-tidy, formatting, and unit tests are mandatory
- Large features require an ADR (architecture decision record)

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

EdgeFEM builds on excellent open-source projects:
- [Eigen](https://eigen.tuxfamily.org/) - Linear algebra
- [Gmsh](https://gmsh.info/) - Mesh generation
- [pybind11](https://github.com/pybind/pybind11) - Python bindings
