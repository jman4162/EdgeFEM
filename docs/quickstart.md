# Quick Start

This guide gets you running your first EdgeFEM simulation in 5 minutes.

## Prerequisites

- EdgeFEM built with Python SDK (see [Installation](installation.md))
- Python path configured: `export PYTHONPATH="build/python:python:$PYTHONPATH"`

## Your First Simulation: Waveguide S-Parameters

Let's compute the S-parameters of a simple WR-90 rectangular waveguide at 10 GHz.

### Using High-Level Design Class

```python
from edgefem.designs import RectWaveguideDesign
import numpy as np

# Create WR-90 waveguide (X-band standard)
wg = RectWaveguideDesign(
    a=22.86e-3,   # Width: 22.86 mm (broadwall)
    b=10.16e-3,   # Height: 10.16 mm (narrowwall)
    length=50e-3  # Length: 50 mm
)

# Print design info
print(f"Cutoff frequency: {wg.cutoff_frequency/1e9:.3f} GHz")

# Generate mesh (Gmsh creates tetrahedral mesh automatically)
wg.generate_mesh(density=10)  # 10 elements per wavelength

# Compute S-parameters at 10 GHz
S = wg.sparams_at_freq(10e9)

# Display results
print(f"\nS-parameters at 10 GHz:")
print(f"  S11: {S[0,0]:.4f} ({20*np.log10(abs(S[0,0])+1e-10):.1f} dB)")
print(f"  S21: {S[1,0]:.4f} ({20*np.log10(abs(S[1,0])):.1f} dB)")
print(f"  |S21| phase: {np.angle(S[1,0])*180/np.pi:.1f} deg")
```

**Expected output:**
```
Cutoff frequency: 6.557 GHz

S-parameters at 10 GHz:
  S11: 0.0012+0.0008j (-56.2 dB)
  S21: 0.9998-0.0156j (-0.0 dB)
  |S21| phase: -0.9 deg
```

### Frequency Sweep

```python
# Sweep from 8 to 12 GHz
freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=21, verbose=True)

# Export to Touchstone format
wg.export_touchstone("waveguide_sweep", freqs, S_list)
print("Saved: waveguide_sweep.s2p")
```

### Plotting Results

```python
import edgefem.plots as vp

# Plot S-parameters vs frequency
vp.plot_sparams_vs_freq(
    freqs, S_list,
    params=['S11', 'S21'],
    title="WR-90 Waveguide S-Parameters",
    save="waveguide_sparams.png"
)
```

## Patch Antenna Example

Design a 2.4 GHz microstrip patch antenna:

```python
from edgefem.designs import PatchAntennaDesign
import numpy as np

# Create patch antenna on FR-4 substrate
patch = PatchAntennaDesign(
    patch_length=29e-3,      # ~lambda/2 at 2.4 GHz
    patch_width=38e-3,       # Width for 50 ohm match
    substrate_height=1.6e-3, # Standard FR-4 thickness
    substrate_eps_r=4.4,     # FR-4 permittivity
)

# Configure probe feed
patch.set_probe_feed(x_offset=0, y_offset=-5e-3)

# Generate mesh
patch.generate_mesh(density=15)

# Compute input impedance at design frequency
Z_in = patch.input_impedance(2.4e9)
print(f"Input impedance: {Z_in.real:.1f} + j{Z_in.imag:.1f} ohms")

# Return loss
rl = patch.return_loss(2.4e9)
print(f"Return loss: {rl:.1f} dB")

# Radiation pattern
pattern = patch.radiation_pattern(2.4e9, n_theta=91, n_phi=73)
D = patch.directivity(2.4e9)
print(f"Directivity: {10*np.log10(D):.1f} dBi")
```

## Unit Cell / Metasurface Example

Analyze a periodic unit cell:

```python
from edgefem.designs import UnitCellDesign

# Create unit cell
cell = UnitCellDesign(
    period_x=5e-3,
    period_y=5e-3,
    substrate_height=1.0e-3,
    substrate_eps_r=3.5,
)

# Add square patch
cell.add_patch(width=4e-3, length=4e-3)

# Generate mesh with periodic boundaries
cell.generate_mesh(density=20)

# Compute reflection coefficient at normal incidence
R, T = cell.reflection_transmission(10e9, theta=0, phi=0, pol='TE')
print(f"Reflection: {abs(R):.4f} ({20*np.log10(abs(R)+1e-10):.1f} dB)")

# Surface impedance
Zs = cell.surface_impedance(10e9)
print(f"Surface impedance: {Zs.real:.1f} + j{Zs.imag:.1f} ohms")
```

## Low-Level API Example

For more control, use the core API directly:

```python
import pyedgefem as em
import numpy as np

# Load mesh from Gmsh file
mesh = em.load_gmsh("examples/rect_waveguide.msh")
print(f"Mesh: {len(mesh.nodes)} nodes, {len(mesh.elements)} elements")

# Build boundary conditions (PEC on tag 1)
bc = em.build_edge_pec(mesh, 1)

# Set up Maxwell parameters
params = em.MaxwellParams()
params.omega = 2 * np.pi * 10e9  # 10 GHz
params.port_abc_scale = 0.5

# Build ports (extract from physical tags 2 and 3)
port_dim = em.RectWaveguidePort()
port_dim.a = 22.86e-3
port_dim.b = 10.16e-3

mode = em.solve_te10_mode(port_dim, 10e9)

surface1 = em.extract_surface_mesh(mesh, 2)
em.populate_te10_field(surface1, port_dim, mode)
port1 = em.build_wave_port(mesh, surface1, mode)

surface2 = em.extract_surface_mesh(mesh, 3)
em.populate_te10_field(surface2, port_dim, mode)
port2 = em.build_wave_port(mesh, surface2, mode)

ports = [port1, port2]

# Compute S-parameters
S = em.calculate_sparams_eigenmode(mesh, params, bc, ports)
print(f"S11 = {S[0,0]:.4f}")
print(f"S21 = {S[1,0]:.4f}")
```

## VTK Export for ParaView

Export fields for 3D visualization in ParaView:

```python
# Export E-field to VTK
wg.export_fields_vtk("waveguide_field.vtu", freq=10e9, active_port=0)
print("Open waveguide_field.vtu in ParaView to visualize E-field")
```

## Next Steps

- **[Waveguide Tutorial](tutorials/waveguide.md)** - In-depth waveguide analysis
- **[Patch Antenna Tutorial](tutorials/patch_antenna.md)** - Antenna design workflow
- **[Unit Cell Tutorial](tutorials/unit_cell.md)** - Metasurface characterization
- **[API Reference](api/core.md)** - Complete function documentation
