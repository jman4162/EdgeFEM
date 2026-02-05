# Installation

VectorEM can be installed from source. Pre-built packages are planned for future releases.

## Prerequisites

### macOS (Apple Silicon & Intel)

```bash
# Install Xcode Command Line Tools
xcode-select --install

# Install dependencies via Homebrew
brew install cmake ninja eigen gmsh
```

### Linux (Ubuntu/Debian)

```bash
sudo apt update
sudo apt install cmake ninja-build libeigen3-dev gmsh python3-dev
```

### Linux (Fedora/RHEL)

```bash
sudo dnf install cmake ninja-build eigen3-devel gmsh python3-devel
```

### Windows

1. Install [Visual Studio 2022](https://visualstudio.microsoft.com/) with C++ workload
2. Install [CMake](https://cmake.org/download/) (3.20+)
3. Install [vcpkg](https://vcpkg.io/) and use it to install Eigen3
4. Install [Gmsh](https://gmsh.info/) and add to PATH

## Building from Source

### Basic Build (C++ only)

```bash
git clone https://github.com/jman4162/VectorEM.git
cd VectorEM

# Configure and build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### With Python SDK

```bash
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DVECTOREM_PYTHON=ON

cmake --build build -j
```

The Python module will be built as `build/python/pyvectorem.cpython-*.so`.

### Verify Installation

```bash
# Run tests
ctest --test-dir build -j

# Test Python bindings
python3 -c "import sys; sys.path.insert(0, 'build/python'); import pyvectorem as em; print('OK')"
```

## Python Environment Setup

For development, we recommend using a virtual environment:

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install Python dependencies
pip install numpy matplotlib

# Optional: for interactive 3D plots
pip install plotly

# Add VectorEM to Python path (development)
export PYTHONPATH="$PWD/build/python:$PWD/python:$PYTHONPATH"
```

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `VECTOREM_BUILD_SCALAR` | ON | Build scalar Helmholtz solver |
| `VECTOREM_BUILD_VECTOR` | ON | Build Maxwell vector solver |
| `VECTOREM_PYTHON` | OFF | Build Python bindings (pybind11) |
| `CMAKE_BUILD_TYPE` | Release | Build type (Release/Debug/RelWithDebInfo) |

Example with all options:

```bash
cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DVECTOREM_BUILD_SCALAR=ON \
    -DVECTOREM_BUILD_VECTOR=ON \
    -DVECTOREM_PYTHON=ON
```

## Troubleshooting

### Gmsh not found

Ensure Gmsh is installed and in your PATH:

```bash
which gmsh  # Should print path
gmsh --version  # Should print version
```

If using a custom Gmsh installation, set the path:

```bash
export PATH="/path/to/gmsh/bin:$PATH"
```

### Eigen not found

If CMake can't find Eigen, specify the path:

```bash
cmake -S . -B build -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 ...
```

### Python module import fails

1. Ensure the module is built: check for `pyvectorem*.so` in `build/python/`
2. Add to Python path: `export PYTHONPATH="$PWD/build/python:$PYTHONPATH"`
3. Check Python version matches the one used during build

### Build fails on Apple Silicon

Ensure you're using native ARM64 tools:

```bash
# Check architecture
uname -m  # Should print "arm64"

# Verify CMake is ARM64
file $(which cmake)  # Should include "arm64"
```

If using Rosetta, rebuild with native tools.

## Next Steps

- [Quick Start](quickstart.md) - Run your first simulation
- [Tutorials](tutorials/waveguide.md) - Step-by-step examples
