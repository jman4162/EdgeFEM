# VectorEM

**3D Finite-Element Electromagnetics Simulator for RF and mmWave**

VectorEM is an open-source full-wave FEM solver specialized for metasurface, metamaterial, and phased array unit cell modeling (100 kHz - 110 GHz). It computes S-parameters, electromagnetic fields, and radiation patterns using Nédélec (edge) elements with industry-standard accuracy.

## Key Features

| Feature | Description |
|---------|-------------|
| **S-parameters** | Wave ports with eigenmode-based extraction (99.4% accuracy) |
| **Radiation Patterns** | 3D far-field via Stratton-Chu integration |
| **Periodic Structures** | Floquet ports for infinite array/metasurface analysis |
| **Open Boundaries** | PML and ABC for antenna/scattering problems |
| **Python SDK** | Full scriptable API with NumPy integration |
| **High-level Design Classes** | `RectWaveguideDesign`, `PatchAntennaDesign`, `UnitCellDesign` |
| **Standard Exports** | Touchstone (.sNp), VTK, CSV |

## Quick Example

```python
from vectorem.designs import RectWaveguideDesign
import numpy as np

# Create WR-90 waveguide (X-band: 8.2-12.4 GHz)
wg = RectWaveguideDesign(
    a=22.86e-3,   # Width: 22.86 mm
    b=10.16e-3,   # Height: 10.16 mm
    length=50e-3  # Length: 50 mm
)

# Generate mesh and compute S-parameters at 10 GHz
wg.generate_mesh(density=10)
S = wg.sparams_at_freq(10e9)

print(f"|S11| = {abs(S[0,0]):.4f}")  # ~0.0 (matched)
print(f"|S21| = {abs(S[1,0]):.4f}")  # ~1.0 (full transmission)
```

## Design Philosophy

VectorEM is designed as a **computation engine** that integrates with your existing workflow:

- **Script-first**: All simulations are fully reproducible Python scripts
- **Standard formats**: Export to Touchstone for RF tools, VTK for ParaView
- **Researcher-friendly**: No GUI lock-in, runs on laptops to HPC clusters
- **Accurate**: Validated against analytical solutions and commercial tools

## Architecture

```
                    ┌─────────────────────────────────┐
                    │      Python SDK (pyvectorem)    │
                    │  ┌───────────┐ ┌─────────────┐  │
                    │  │  Design   │ │   Plots     │  │
                    │  │  Classes  │ │   Module    │  │
                    │  └───────────┘ └─────────────┘  │
                    └─────────────┬───────────────────┘
                                  │
┌─────────────────────────────────┼─────────────────────────────────┐
│                     C++ Core    │                                 │
│  ┌──────────┐ ┌──────────┐ ┌───┴────┐ ┌──────────┐ ┌──────────┐  │
│  │   Mesh   │ │  Maxwell │ │ Ports  │ │  Solver  │ │  Post-   │  │
│  │  (Gmsh)  │ │ Assembly │ │(Wave/  │ │(BiCGSTAB)│ │  proc    │  │
│  │          │ │ (PML/ABC)│ │Floquet)│ │          │ │(NTF,VTK) │  │
│  └──────────┘ └──────────┘ └────────┘ └──────────┘ └──────────┘  │
└───────────────────────────────────────────────────────────────────┘
```

## Getting Started

1. **[Installation](installation.md)** - Set up VectorEM on your system
2. **[Quick Start](quickstart.md)** - Run your first simulation in 5 minutes
3. **[Tutorials](tutorials/waveguide.md)** - Step-by-step guides for common tasks

## Ecosystem Integration

VectorEM is the full-wave FEM engine in a multi-package RF modeling ecosystem:

- **[Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model)**: Uses VectorEM element patterns and coupling for array synthesis
- **Metasurface Package** (planned): Unit cell optimization and homogenization

## License

VectorEM is released under the MIT License. See `LICENSE` for details.
