# EdgeFEM

**3D Finite-Element Electromagnetics Simulator for RF and mmWave**

EdgeFEM is an open-source full-wave frequency-domain FEM solver for RF and microwave structures: hollow waveguides, probe-fed patch antennas, and periodic unit cells. It computes S-parameters, electromagnetic fields, and radiation patterns using lowest-order Nédélec (edge) tetrahedral elements. Accuracy is validated against analytical benchmarks over roughly 5-26 GHz (see [Validation](validation.md)).

## Key Features

| Feature | Description |
|---------|-------------|
| **S-parameters** | Wave ports with 2D discrete port modes and assembled surface-mass ABC; lumped ports for antenna feeds |
| **Radiation Patterns** | 3D far-field via Stratton-Chu integration |
| **Periodic Structures** | Phase-shifted (Bloch) periodic BCs along one lattice axis, for small unit cells |
| **Open Boundaries** | First-order ABC and a graded absorbing layer |
| **Python SDK** | Full scriptable API with NumPy integration |
| **High-level Design Classes** | `RectWaveguideDesign`, `PatchAntennaDesign`, `UnitCellDesign` |
| **Standard Exports** | Touchstone (.sNp), VTK, CSV |

## Quick Example

```python
from edgefem.designs import RectWaveguideDesign
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

EdgeFEM is designed as a **computation engine** that integrates with your existing workflow:

- **Script-first**: All simulations are fully reproducible Python scripts
- **Standard formats**: Export to Touchstone for RF tools, VTK for ParaView
- **Researcher-friendly**: No GUI lock-in; single-node, script-driven workflow
- **Accurate**: Validated against analytical solutions (waveguide dispersion, cavity eigenmodes, Fresnel)

## Architecture

```
                    ┌─────────────────────────────────┐
                    │      Python SDK (pyedgefem)    │
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
│  │          │ │ (ABC/abs)│ │Lumped) │ │          │ │(NTF,VTK) │  │
│  └──────────┘ └──────────┘ └────────┘ └──────────┘ └──────────┘  │
└───────────────────────────────────────────────────────────────────┘
```

## Getting Started

1. **[Installation](installation.md)** - Set up EdgeFEM on your system
2. **[Quick Start](quickstart.md)** - Run your first simulation in 5 minutes
3. **[Tutorials](tutorials/waveguide.md)** - Step-by-step guides for common tasks

## Ecosystem Integration

EdgeFEM is the full-wave FEM engine in a multi-package RF modeling ecosystem:

- **[Phased-Array-Antenna-Model](https://github.com/jman4162/Phased-Array-Antenna-Model)**: Uses EdgeFEM element patterns and coupling for array synthesis
- **Metasurface Package** (planned): Unit cell optimization and homogenization

## License

EdgeFEM is released under the MIT License. See `LICENSE` for details.
