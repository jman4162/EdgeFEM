# VectorEM

**3D Finite-Element Electromagnetics Simulator for RF and mmWave**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)

VectorEM is an open-source full-wave FEM solver specialized for metasurface, metamaterial, and phased array unit cell modeling (100 kHz - 110 GHz). It computes S-parameters, electromagnetic fields, and radiation patterns using Nédélec (edge) elements with industry-standard accuracy.

## Features

| Feature | Description |
|---------|-------------|
| **S-parameters** | Wave ports with eigenmode-based extraction (99.4% accuracy) |
| **Radiation Patterns** | 3D far-field via Stratton-Chu integration |
| **Periodic Structures** | Floquet ports for infinite array/metasurface analysis |
| **Open Boundaries** | PML and ABC for antenna/scattering problems |
| **Python SDK** | Full scriptable API with NumPy integration |
| **Design Classes** | `RectWaveguideDesign`, `PatchAntennaDesign`, `UnitCellDesign` |
| **Plotting Module** | 14 publication-quality visualization functions |
| **Standard Exports** | Touchstone (.sNp), VTK, CSV |

## Quick Start

### Installation

```bash
# macOS
brew install cmake ninja eigen gmsh

# Clone and build
git clone https://github.com/jman4162/VectorEM.git
cd VectorEM
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DVECTOREM_PYTHON=ON
cmake --build build -j

# Set Python path
export PYTHONPATH="$PWD/build/python:$PWD/python:$PYTHONPATH"
```

### 10-Line Example: Waveguide S-Parameters

```python
from vectorem.designs import RectWaveguideDesign
import numpy as np

# Create WR-90 waveguide (X-band)
wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
wg.generate_mesh(density=10)

# Compute S-parameters at 10 GHz
S = wg.sparams_at_freq(10e9)
print(f"|S11| = {abs(S[0,0]):.4f}")  # ~0.001 (matched)
print(f"|S21| = {abs(S[1,0]):.4f}")  # ~0.994 (99.4% transmission)
```

### Patch Antenna Example

```python
from vectorem.designs import PatchAntennaDesign

# 2.4 GHz WiFi patch on FR-4
patch = PatchAntennaDesign(
    patch_length=29e-3, patch_width=38e-3,
    substrate_height=1.6e-3, substrate_eps_r=4.4
)
patch.set_probe_feed(y_offset=-5e-3)
patch.generate_mesh(density=15)

# Analyze
Z_in = patch.input_impedance(2.4e9)
pattern = patch.radiation_pattern(2.4e9)
D = patch.directivity(2.4e9)
print(f"Directivity: {10*np.log10(D):.1f} dBi")
```

### Unit Cell / Metasurface Example

```python
from vectorem.designs import UnitCellDesign

# Periodic unit cell
cell = UnitCellDesign(period_x=5e-3, period_y=5e-3,
                       substrate_height=0.5e-3, substrate_eps_r=3.5)
cell.add_patch(width=4e-3, length=4e-3)
cell.generate_mesh(density=20)

# Reflection coefficient at normal incidence
R, T = cell.reflection_transmission(10e9, theta=0, phi=0)
Zs = cell.surface_impedance(10e9)
```

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
│                      Python SDK (pyvectorem)                     │
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
│  │          │ │(PML/ABC) │ │Floquet)│ │          │ │  VTK)   │  │
│  └──────────┘ └──────────┘ └────────┘ └──────────┘ └─────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

## Validation

VectorEM achieves research-grade accuracy:

| Metric | Target | Achieved |
|--------|--------|----------|
| \|S21\| transmission | > 99% | 99.4% |
| Phase accuracy | < 5° | 1-2° |
| Passivity | \|S11\|² + \|S21\|² ≤ 1 | PASS |
| Reciprocity | S12 = S21 | PASS |

See [Validation Report](docs/validation.md) for full benchmark results.

## Ecosystem Integration

VectorEM is the full-wave FEM engine in a multi-package RF modeling ecosystem:

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
│                          VectorEM                                │
│  Full-wave 3D FEM providing: S-params, patterns, coupling       │
└─────────────────────────────────────────────────────────────────┘
```

- **[Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model)**: Uses VectorEM element patterns and coupling for array synthesis

## Build Options

```bash
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DVECTOREM_PYTHON=ON         # Python bindings
```

| Option | Default | Description |
|--------|---------|-------------|
| `VECTOREM_BUILD_SCALAR` | ON | Scalar Helmholtz solver |
| `VECTOREM_BUILD_VECTOR` | ON | Maxwell vector solver |
| `VECTOREM_PYTHON` | OFF | Python bindings (pybind11) |

## Run Tests

```bash
ctest --test-dir build -j        # All tests
ctest --test-dir build -L smoke  # Quick smoke tests
ctest --test-dir build -L pml    # PML absorption tests
```

## Plotting Module

```python
import vectorem.plots as vp

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

If you use VectorEM in your research, please cite:

```bibtex
@software{hodge2026vectorem,
  author = {Hodge, John},
  title = {{VectorEM}: A 3D Finite-Element Electromagnetics Simulator for RF and mmWave},
  year = {2026},
  url = {https://github.com/jman4162/VectorEM},
  version = {1.0.0}
}
```

## Contributing

- Use feature branches; PRs require green CI and updated docs
- clang-tidy, formatting, and unit tests are mandatory
- Large features require an ADR (architecture decision record)

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

VectorEM builds on excellent open-source projects:
- [Eigen](https://eigen.tuxfamily.org/) - Linear algebra
- [Gmsh](https://gmsh.info/) - Mesh generation
- [pybind11](https://github.com/pybind/pybind11) - Python bindings
