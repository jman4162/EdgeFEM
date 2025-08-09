# VectorEM

VectorEM is a 3D finite‑element electromagnetics simulator for RF and mmWave. It computes S‑parameters, fields, and antenna radiation patterns using Nédélec (edge) elements, driven‑modal/terminal ports, and PML for open boundaries. It scales from laptops to HPC clusters and exposes a Python SDK for scripting.

## Features (v1)
- Frequency‑domain Maxwell FEM with 1st–3rd order edge elements
- Wave and lumped ports; S‑parameters with Touchstone export
- PML radiation boundaries; ABC fallback
- Adaptive h/p refinement with goal‑oriented estimators
- Far‑field patterns via Huygens surface / Stratton–Chu
- Iterative (GMRES+AMG) and optional direct solvers; MPI support
- Python API and headless CLI

## Repo Layout
```

include/vectorem/   # public headers
src/                # implementation
examples/           # meshes & scripts
python/             # pybind11 bindings
cmake/              # cmake helpers

````

## Quick Start

### Prereqs
- C++20 toolchain, CMake ≥ 3.20
- Eigen ≥ 3.4
- (Optional) Python 3.11, pybind11 for SDK
- macOS: Xcode Command Line Tools (`xcode-select --install`), Homebrew (`brew install cmake ninja eigen`)

### Build
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
````

### macOS (Apple Silicon) quick setup

```bash
# One-time
xcode-select --install || true
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" || true
brew update
brew install cmake ninja eigen

# Build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

> Note: v1 uses CPU on macOS; NVIDIA GPU acceleration is Linux-only. Metal/MPS backend is planned for v1.x.

### Run the Scalar Demo (Milestone 0)

```bash
./build/vectorem_scalar_demo examples/cube_cavity.msh 1
```

### Run Tests

```bash
ctest --test-dir build -j
# PML absorption tests (skipped by default)
ctest --test-dir build -L pml
```

## Python SDK (optional)

```bash
cmake -S . -B build -DVECTOREM_PYTHON=ON
cmake --build build -j
python -c "import pyvectorem as em; print('ok')"
```

## Design Highlights

* Nédélec edge elements; auxiliary‑space AMG for robust preconditioning
* Ports via 2D eigenmodes; rigorous S‑param normalization
* PML with automatic thickness/grading; diagnostics for reflections

## Validation

See `docs/validation.md` for the golden benchmark suite and tolerances.

## Contributing

* Use feature branches; PRs require green CI and updated docs
* clang‑tidy, formatting, and unit tests are mandatory
* Large features require an ADR (architecture decision record)

## License

TBD. See `LICENSE`.
