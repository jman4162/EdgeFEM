# Wave Port Workflow Guide

This guide explains how to set up and extract S-parameters using modal wave ports in VectorEM.

## Overview

Wave ports provide accurate S-parameter extraction for waveguide-based structures. The workflow involves:

1. Creating a mesh with port surfaces defined as physical groups
2. Extracting the 2D port surface mesh
3. Solving the 2D eigenmode problem on the port cross-section
4. Projecting the modal field onto 3D edge elements
5. Assembling the Maxwell system with port boundary conditions
6. Solving and extracting S-parameters

## Prerequisites

- VectorEM built with `VECTOREM_BUILD_VECTOR=ON`
- Gmsh for mesh generation
- Python bindings (optional): build with `VECTOREM_PYTHON=ON`

## Step 1: Create Geometry with Port Surfaces

In your Gmsh geometry file, define port surfaces as physical groups:

```geo
// rect_waveguide.geo
// WR-90 waveguide: a=22.86mm, b=10.16mm, length=50mm

a = 0.02286;
b = 0.01016;
L = 0.05;

// Create waveguide box
Point(1) = {0, 0, 0};
Point(2) = {a, 0, 0};
Point(3) = {a, b, 0};
Point(4) = {0, b, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Extrude to create volume
Extrude {0, 0, L} { Surface{1}; }

// Physical groups
Physical Surface("PEC", 1) = {13, 17, 21, 25};  // Walls
Physical Surface("Port1", 2) = {1};              // Input port (z=0)
Physical Surface("Port2", 3) = {26};             // Output port (z=L)
Physical Volume("Air", 100) = {1};
```

Generate the mesh:
```bash
gmsh rect_waveguide.geo -3 -o rect_waveguide.msh
```

## Step 2: Load Mesh and Define Boundaries

### C++ API

```cpp
#include "vectorem/mesh.hpp"
#include "vectorem/bc.hpp"
#include "vectorem/maxwell.hpp"
#include "vectorem/ports/wave_port.hpp"
#include "vectorem/ports/port_eigensolve.hpp"

using namespace vectorem;

// Load mesh
Mesh mesh = load_gmsh_v2("rect_waveguide.msh");

// Apply PEC boundary conditions to waveguide walls
BC bc = build_edge_pec(mesh, 1);  // Physical tag 1 = PEC
```

### Python API

```python
import pyvectorem as em

mesh = em.load_gmsh("rect_waveguide.msh")
bc = em.build_edge_pec(mesh, 1)
```

## Step 3: Extract Port Surface Meshes

The port surface mesh is extracted from the volume mesh using the port's physical tag.

### C++ API

```cpp
int port1_tag = 2;
int port2_tag = 3;

PortSurfaceMesh port1_surface = extract_surface_mesh(mesh, port1_tag);
PortSurfaceMesh port2_surface = extract_surface_mesh(mesh, port2_tag);
```

### Python API

```python
port1_surface = em.extract_surface_mesh(mesh, 2)
port2_surface = em.extract_surface_mesh(mesh, 3)
```

## Step 4: Solve Port Eigenmodes

Compute the 2D modal field distribution on each port cross-section.

### C++ API

```cpp
double freq = 10e9;  // Operating frequency
double omega = 2 * M_PI * freq;
std::complex<double> eps_r = 1.0;
std::complex<double> mu_r = 1.0;

// Solve for TE modes (dominant mode in rectangular waveguide)
auto modes1 = solve_port_eigens(port1_surface.mesh, 1, omega, eps_r, mu_r,
                                 ModePolarization::TE);
auto modes2 = solve_port_eigens(port2_surface.mesh, 1, omega, eps_r, mu_r,
                                 ModePolarization::TE);

PortMode mode1 = modes1[0];  // TE10 mode
PortMode mode2 = modes2[0];
```

### Python API

```python
import numpy as np

freq = 10e9
omega = 2 * np.pi * freq
eps_r = 1.0 + 0j
mu_r = 1.0 + 0j

modes1 = em.solve_port_eigens(port1_surface.mesh, 1, omega, eps_r, mu_r,
                               em.ModePolarization.TE)
modes2 = em.solve_port_eigens(port2_surface.mesh, 1, omega, eps_r, mu_r,
                               em.ModePolarization.TE)

mode1 = modes1[0]
mode2 = modes2[0]
```

### Modal Properties

The `PortMode` structure contains:
- `fc`: Cutoff frequency (Hz)
- `kc`: Cutoff wavenumber (rad/m)
- `beta`: Propagation constant
- `Z0`: Modal impedance (Ohms)
- `field`: 2D modal field samples

For WR-90 waveguide, the TE10 mode has:
- Cutoff: fc = c/(2a) ≈ 6.56 GHz
- Modal impedance: Z0 = η₀/√(1 - (fc/f)²)

## Step 5: Build Wave Ports

Project the 2D modal field onto the 3D edge elements at the port surface.

### C++ API

```cpp
WavePort wave_port1 = build_wave_port(mesh, port1_surface, mode1);
WavePort wave_port2 = build_wave_port(mesh, port2_surface, mode2);

std::vector<WavePort> ports = {wave_port1, wave_port2};
```

### Python API

```python
wave_port1 = em.build_wave_port(mesh, port1_surface, mode1)
wave_port2 = em.build_wave_port(mesh, port2_surface, mode2)
ports = [wave_port1, wave_port2]
```

## Step 6: Assemble Maxwell System and Calculate S-Parameters

### C++ API

```cpp
MaxwellParams params;
params.omega = omega;
params.eps_r = eps_r;
params.mu_r = mu_r;

// Calculate full S-parameter matrix
Eigen::MatrixXcd S = calculate_sparams(mesh, params, bc, ports);

// S is now a 2x2 complex matrix
std::cout << "S11 = " << S(0,0) << std::endl;
std::cout << "S21 = " << S(1,0) << std::endl;
```

### Python API

```python
params = em.MaxwellParams()
params.omega = omega
params.eps_r = eps_r
params.mu_r = mu_r

S = em.calculate_sparams(mesh, params, bc, ports)

print(f"S11 = {S[0,0]:.4f} ({20*np.log10(abs(S[0,0])):.2f} dB)")
print(f"S21 = {S[1,0]:.4f} ({20*np.log10(abs(S[1,0])):.2f} dB)")
```

## Step 7: Frequency Sweep

For broadband characterization, repeat the calculation at multiple frequencies.

```python
freqs = np.linspace(8e9, 12e9, 51)
s11_data = []
s21_data = []

for f in freqs:
    params.omega = 2 * np.pi * f

    # Re-solve modes at each frequency for accurate dispersion
    modes1 = em.solve_port_eigens(port1_surface.mesh, 1, params.omega,
                                   eps_r, mu_r, em.ModePolarization.TE)
    modes2 = em.solve_port_eigens(port2_surface.mesh, 1, params.omega,
                                   eps_r, mu_r, em.ModePolarization.TE)

    wp1 = em.build_wave_port(mesh, port1_surface, modes1[0])
    wp2 = em.build_wave_port(mesh, port2_surface, modes2[0])

    S = em.calculate_sparams(mesh, params, bc, [wp1, wp2])
    s11_data.append(S[0,0])
    s21_data.append(S[1,0])
```

## Step 8: Export to Touchstone

```python
sparams_list = []
for i in range(len(freqs)):
    s = em.SParams2()
    s.s11 = s11_data[i]
    s.s21 = s21_data[i]
    s.s12 = s21_data[i]  # Reciprocal
    s.s22 = s11_data[i]  # Symmetric
    sparams_list.append(s)

em.write_touchstone("waveguide.s2p", list(freqs), sparams_list)
```

## Validation Checks

### Passivity Check

For a passive network: |S11|² + |S21|² ≤ 1

```python
power_sum = np.abs(S[0,0])**2 + np.abs(S[1,0])**2
assert power_sum <= 1.0 + 1e-6, "Network is non-passive"
```

### Reciprocity Check

For a reciprocal network: S12 = S21

```python
assert np.isclose(S[0,1], S[1,0], rtol=1e-4), "Network is non-reciprocal"
```

## Troubleshooting

### Mode Not Found

If `solve_port_eigens` returns an empty list:
- Verify the operating frequency is above the mode's cutoff frequency
- Check that the port surface mesh has the correct physical tag
- Ensure the port cross-section is properly meshed

### Poor S-Parameter Accuracy

- Increase mesh density, especially at port surfaces
- Verify PEC boundary conditions are applied to all walls
- Check that the port surfaces are planar and properly aligned

### Solver Convergence Issues

Use the extended solver options:
```python
opts = em.SolveOptions()
opts.tolerance = 1e-12
opts.max_iterations = 20000
result = em.solve_linear(A, b, opts)
if not result.converged:
    print(f"Warning: {result.error_message}")
```

## Example Files

- `examples/rect_waveguide.geo` - Gmsh geometry for WR-90 waveguide
- `examples/rect_waveguide.msh` - Pre-generated mesh
- `examples/waveport_demo.py` - Complete Python example
- `examples/waveguide_sparams.s2p` - Reference S-parameter data
