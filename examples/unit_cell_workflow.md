# Unit Cell Simulation Workflow

This guide explains how to simulate phased array unit cells using periodic boundary conditions (Floquet analysis) in EdgeFEM.

## Overview

Periodic boundary conditions enable efficient simulation of infinite arrays by modeling a single unit cell with phase-shifted boundary conditions. This approach:

- Reduces computational cost from O(N²) elements to O(1)
- Captures mutual coupling effects through Floquet modes
- Enables scan angle analysis via phase shift variation

## Prerequisites

- EdgeFEM built with `EDGEFEM_BUILD_VECTOR=ON`
- Gmsh for mesh generation
- Python bindings (recommended): `EDGEFEM_PYTHON=ON`

## Step 1: Create Unit Cell Geometry

Design your unit cell in Gmsh with clearly defined periodic surface pairs.

```geo
// unit_cell.geo
cell_size = 0.015;  // 15mm (~λ/2 at 10 GHz)

// Define unit cell boundaries
// ... geometry definition ...

// Physical groups for periodic BC
Physical Surface("MasterX", 4) = {left_face};   // x = 0
Physical Surface("SlaveX", 5) = {right_face};   // x = cell_size
Physical Surface("MasterY", 6) = {front_face};  // y = 0
Physical Surface("SlaveY", 7) = {back_face};    // y = cell_size

// Floquet port (radiation boundary)
Physical Surface("FloquetPort", 3) = {top_face};

// PEC boundaries
Physical Surface("Ground", 1) = {bottom_face};
Physical Surface("Patch", 2) = {patch_surface};
```

Generate the mesh:
```bash
gmsh unit_cell.geo -3 -o unit_cell.msh
```

## Step 2: Set Up Periodic Boundary Conditions

### C++ API

```cpp
#include "edgefem/periodic.hpp"
#include "edgefem/maxwell.hpp"

using namespace edgefem;

// Load mesh
Mesh mesh = load_gmsh_v2("unit_cell.msh");

// Define period vectors
Eigen::Vector3d period_x(cell_size, 0, 0);
Eigen::Vector3d period_y(0, cell_size, 0);

// Build periodic pairs
PeriodicBC pbc_x = build_periodic_pairs(mesh, 4, 5, period_x);
PeriodicBC pbc_y = build_periodic_pairs(mesh, 6, 7, period_y);

// Validate
assert(validate_periodic_bc(mesh, pbc_x));
assert(validate_periodic_bc(mesh, pbc_y));
```

### Python API

```python
import pyedgefem as em
import numpy as np

mesh = em.load_gmsh("unit_cell.msh")

cell_size = 0.015  # 15mm
period_x = np.array([cell_size, 0, 0])
period_y = np.array([0, cell_size, 0])

pbc_x = em.build_periodic_pairs(mesh, 4, 5, period_x)
pbc_y = em.build_periodic_pairs(mesh, 6, 7, period_y)
```

## Step 3: Configure Floquet Phase Shift

The phase shift determines the scan angle of the array:

```
phase_x = exp(j * kx * Lx)
phase_y = exp(j * ky * Ly)

where:
  kx = k0 * sin(θ) * cos(φ)
  ky = k0 * sin(θ) * sin(φ)
  k0 = 2π/λ (free-space wavenumber)
  θ = scan angle from broadside
  φ = azimuth angle
```

### C++ API

```cpp
double freq = 10e9;
double k0 = 2 * M_PI * freq / 299792458.0;
double theta = 30.0 * M_PI / 180.0;  // 30° scan
double phi = 0.0;  // E-plane

// Compute phase shifts
auto phase_x = floquet_phase_from_angle(period_x, theta, phi, k0);
auto phase_y = floquet_phase_from_angle(period_y, theta, phi, k0);

// Apply to periodic BC
set_floquet_phase(pbc_x, Eigen::Vector2d(k0 * sin(theta) * cos(phi), 0));
```

### Python API

```python
import numpy as np

c0 = 299792458.0
freq = 10e9
k0 = 2 * np.pi * freq / c0
theta = np.radians(30)  # 30° scan
phi = 0  # E-plane

kx = k0 * np.sin(theta) * np.cos(phi)
ky = k0 * np.sin(theta) * np.sin(phi)

em.set_floquet_phase(pbc_x, np.array([kx, ky]))
```

## Step 4: Assemble and Solve with Periodic BC

### C++ API

```cpp
BC bc = build_edge_pec(mesh, 1);  // Ground plane
// Also apply PEC to patch: merge into bc.dirichlet_edges

MaxwellParams params;
params.omega = 2 * M_PI * freq;
params.eps_r = 1.0;
params.mu_r = 1.0;

// Combined periodic BC (merge x and y pairs)
PeriodicBC pbc_combined;
pbc_combined.pairs.insert(pbc_combined.pairs.end(),
                          pbc_x.pairs.begin(), pbc_x.pairs.end());
pbc_combined.pairs.insert(pbc_combined.pairs.end(),
                          pbc_y.pairs.begin(), pbc_y.pairs.end());
pbc_combined.phase_shift = phase_x * phase_y;  // Combined phase

// Set up Floquet port
std::vector<WavePort> ports = { floquet_port };

// Assemble with periodic BC
auto asmbl = assemble_maxwell_periodic(mesh, params, bc, pbc_combined, ports, 0);

// Solve
auto result = solve_linear(asmbl.A, asmbl.b, {});

// Calculate S-parameters
auto S = calculate_sparams_periodic(mesh, params, bc, pbc_combined, ports);
```

### Python API

```python
bc = em.build_edge_pec(mesh, 1)

params = em.MaxwellParams()
params.omega = 2 * np.pi * freq
params.eps_r = 1.0 + 0j
params.mu_r = 1.0 + 0j

# Assemble with periodic BC
asm = em.assemble_maxwell_periodic(mesh, params, bc, pbc_combined, ports, 0)

# Solve
result = em.solve_linear(asm.A, asm.b)
print(f"Converged: {result.converged}, Iterations: {result.iters}")

# Calculate S-parameters
S = em.calculate_sparams_periodic(mesh, params, bc, pbc_combined, ports)
print(f"S11 = {S[0,0]:.4f} ({20*np.log10(abs(S[0,0])):.2f} dB)")
```

## Step 5: Scan Angle Sweep

To characterize active impedance vs. scan angle:

```python
scan_angles = np.linspace(0, 60, 31)  # 0° to 60° in 2° steps
s11_vs_theta = []

for theta_deg in scan_angles:
    theta = np.radians(theta_deg)

    # Update phase shifts
    kx = k0 * np.sin(theta)
    ky = 0  # E-plane scan

    em.set_floquet_phase(pbc_x, np.array([kx, ky]))
    em.set_floquet_phase(pbc_y, np.array([kx, ky]))

    # Recombine PBC
    pbc_combined.phase_shift = np.exp(1j * kx * cell_size)

    # Calculate S-parameters
    S = em.calculate_sparams_periodic(mesh, params, bc, pbc_combined, ports)
    s11_vs_theta.append(S[0,0])

    print(f"θ = {theta_deg:5.1f}°: S11 = {20*np.log10(abs(S[0,0])):+6.2f} dB")
```

## Step 6: Active Impedance Calculation

The active impedance (scan impedance) is related to the active reflection coefficient:

```python
# From S-parameter at scan angle
s11_active = s11_vs_theta[-1]  # At 60° scan

# Active VSWR
vswr = (1 + abs(s11_active)) / (1 - abs(s11_active))
print(f"Active VSWR at 60°: {vswr:.2f}")

# Active impedance
z0 = 50.0
z_active = z0 * (1 + s11_active) / (1 - s11_active)
print(f"Active impedance: {z_active.real:.1f} + j{z_active.imag:.1f} Ω")
```

## Step 7: Export Results

```python
# Export S-parameters vs. frequency at multiple scan angles
coupling_data = []
for theta_deg in [0, 15, 30, 45, 60]:
    # ... compute S-matrix at this angle ...
    cm = em.compute_coupling(S, freq)
    coupling_data.append(cm)

em.write_coupling_json("unit_cell_coupling.json", coupling_data)

# Export Touchstone
S_matrices = [...]  # Collect S-matrices at all frequencies
em.write_touchstone_nport("unit_cell.s1p", list(freqs), S_matrices)
```

## Physical Interpretation

### Active Element Pattern

The embedded element pattern accounts for mutual coupling:
- At broadside (θ=0), elements see nominal impedance
- As scan angle increases, active impedance changes
- Scan blindness can occur when active VSWR → ∞

### Grating Lobes

Grating lobes appear when:
```
d/λ > 1/(1 + sin(θ_max))
```

For λ/2 spacing (d = λ/2):
- First grating lobe at θ = 90° (horizon)
- Safe to scan to ~60° without grating lobes

### Surface Waves

On dielectric substrates, surface waves can cause:
- Scan blindness at specific angles
- Reduced gain at wide angles
- Increased mutual coupling

## Troubleshooting

### No Edge Pairs Found

```
Error: No edges found on master surface with tag 4
```

- Verify physical group tags in Gmsh geometry
- Check that surface triangles have correct physical tags
- Ensure mesh was generated with `-3` (3D) option

### Phase Shift Issues

If S-parameters don't vary smoothly with scan angle:
- Check period vector direction matches surface normals
- Verify edge orientations are consistent
- Ensure phase is applied to correct periodic pair

### Solver Convergence

For unit cells with high-Q resonances:
```python
opts = em.SolveOptions()
opts.tolerance = 1e-12
opts.max_iterations = 20000
opts.use_bicgstab = True
```

## Example Files

- `examples/unit_cell.geo` - Patch element unit cell geometry
- `examples/unit_cell_demo.py` - Complete Python workflow
- `examples/waveport_demo.py` - Wave port reference

## References

1. Munk, B. A., "Frequency Selective Surfaces: Theory and Design"
2. Pozar, D. M., "Microwave Engineering" (Chapter 10)
3. Balanis, C. A., "Antenna Theory: Analysis and Design" (Chapter 6)
