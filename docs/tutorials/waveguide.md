# Tutorial: Rectangular Waveguide S-Parameters

This tutorial covers computing S-parameters of rectangular waveguides using EdgeFEM. You'll learn to:

- Set up waveguide geometry and mesh
- Compute S-parameters at single and multiple frequencies
- Validate results against analytical solutions
- Export data for use in other tools

## Background

Rectangular waveguides are fundamental building blocks in microwave systems. They support TE and TM modes with well-known analytical solutions, making them ideal for validating FEM simulation accuracy.

**Key equations:**

- TE10 cutoff frequency: $f_c = \frac{c}{2a}$ where $a$ is the waveguide width
- Propagation constant: $\beta = \frac{2\pi}{c}\sqrt{f^2 - f_c^2}$
- For a matched waveguide of length $L$: $S_{11} \approx 0$, $S_{21} = e^{-j\beta L}$

## Step 1: Create Waveguide Design

```python
from edgefem.designs import RectWaveguideDesign
import numpy as np

# WR-90 waveguide (X-band: 8.2-12.4 GHz)
a = 22.86e-3  # Width: 22.86 mm
b = 10.16e-3  # Height: 10.16 mm
L = 50e-3     # Length: 50 mm

wg = RectWaveguideDesign(a=a, b=b, length=L)

# Check cutoff frequency
fc = wg.cutoff_frequency
print(f"TE10 cutoff frequency: {fc/1e9:.3f} GHz")
# Expected: ~6.56 GHz
```

## Step 2: Generate Mesh

EdgeFEM uses Gmsh to generate tetrahedral meshes automatically:

```python
# Generate mesh with 10 elements per wavelength
wg.generate_mesh(density=10)

# The mesh is stored internally
print(f"Mesh generated with {len(wg.mesh.elements)} elements")
```

**Mesh density guidelines:**

| Application | Density | Notes |
|-------------|---------|-------|
| Quick check | 5-8 | Fast, ~5% accuracy |
| Standard analysis | 10-15 | Good balance |
| High accuracy | 15-25 | Slower, <1% error |
| Validation | 25+ | Reference results |

To save the mesh for inspection in Gmsh:

```python
wg.generate_mesh(density=10, output_path="waveguide_mesh")
# Creates: waveguide_mesh.geo and waveguide_mesh.msh
```

## Step 3: Single-Frequency Analysis

Compute S-parameters at 10 GHz:

```python
freq = 10e9  # 10 GHz

S = wg.sparams_at_freq(freq)

# Display results
print(f"\nS-parameters at {freq/1e9:.1f} GHz:")
print(f"  S11: {S[0,0]:.6f}")
print(f"  S21: {S[1,0]:.6f}")
print(f"  S12: {S[0,1]:.6f}")
print(f"  S22: {S[1,1]:.6f}")

# Magnitude in dB
print(f"\nMagnitude:")
print(f"  |S11| = {20*np.log10(abs(S[0,0])+1e-10):.1f} dB")
print(f"  |S21| = {20*np.log10(abs(S[1,0])):.2f} dB")
```

**Expected results** for matched 50mm waveguide at 10 GHz:
- |S11| ≈ -50 dB (excellent match)
- |S21| ≈ 0 dB (full transmission)

## Step 4: Compare with Analytical Solution

```python
# Analytical propagation constant
c0 = 299792458.0
beta_analytical = (2 * np.pi / c0) * np.sqrt(freq**2 - fc**2)

# Analytical S21 (lossless matched waveguide)
S21_analytical = np.exp(-1j * beta_analytical * L)

# Compare
S21_sim = S[1,0]
print(f"\nComparison at {freq/1e9:.1f} GHz:")
print(f"  Analytical S21: {S21_analytical:.6f}")
print(f"  Simulated S21:  {S21_sim:.6f}")

# Phase comparison
phase_ana = np.angle(S21_analytical) * 180 / np.pi
phase_sim = np.angle(S21_sim) * 180 / np.pi
print(f"\n  Analytical phase: {phase_ana:.2f} deg")
print(f"  Simulated phase:  {phase_sim:.2f} deg")
print(f"  Phase error:      {abs(phase_ana - phase_sim):.2f} deg")

# Magnitude comparison
mag_error = abs(abs(S21_sim) - 1.0) * 100
print(f"\n  |S21| error: {mag_error:.2f}%")
```

EdgeFEM achieves <1% magnitude error and <2° phase error for typical meshes.

## Step 5: Frequency Sweep

Compute S-parameters across X-band:

```python
# Sweep from 8 to 12 GHz
freqs, S_list = wg.frequency_sweep(
    f_start=8e9,
    f_stop=12e9,
    n_points=41,
    verbose=True
)

# S_list is a list of 2x2 matrices, one per frequency
print(f"\nComputed {len(freqs)} frequency points")
```

## Step 6: Visualize Results

```python
import matplotlib.pyplot as plt
import edgefem.plots as vp

# Using the plots module
vp.plot_sparams_vs_freq(
    freqs, S_list,
    params=['S11', 'S21'],
    show_phase=True,
    title="WR-90 Waveguide S-Parameters",
    save="waveguide_sparams.png"
)
```

Or create custom plots:

```python
# Extract S21 vs frequency
S21_mag = np.array([abs(S[1,0]) for S in S_list])
S21_phase = np.array([np.angle(S[1,0]) * 180/np.pi for S in S_list])

# Analytical comparison
beta = (2 * np.pi / c0) * np.sqrt(freqs**2 - fc**2)
S21_ana_phase = -beta * L * 180 / np.pi

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Magnitude
axes[0].plot(freqs/1e9, 20*np.log10(S21_mag), 'b-', label='Simulated')
axes[0].axhline(y=0, color='r', linestyle='--', label='Ideal')
axes[0].set_xlabel('Frequency (GHz)')
axes[0].set_ylabel('|S21| (dB)')
axes[0].set_title('WR-90 Waveguide Transmission')
axes[0].legend()
axes[0].grid(True)
axes[0].set_ylim(-1, 0.1)

# Phase
axes[1].plot(freqs/1e9, S21_phase, 'b-', label='Simulated')
axes[1].plot(freqs/1e9, S21_ana_phase, 'r--', label='Analytical')
axes[1].set_xlabel('Frequency (GHz)')
axes[1].set_ylabel('Phase(S21) (deg)')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig('waveguide_validation.png', dpi=150)
plt.show()
```

## Step 7: Export Results

### Touchstone Format

Export for use in circuit simulators:

```python
# Export to Touchstone (.s2p)
wg.export_touchstone("wr90_sweep", freqs, S_list)
print("Saved: wr90_sweep.s2p")
```

The Touchstone file can be imported into:
- AWR Microwave Office
- Keysight ADS
- ANSYS HFSS
- CST Studio
- Any SPICE simulator with S-parameter support

### VTK Export for ParaView

Visualize the E-field distribution:

```python
# Export field at 10 GHz
wg.export_fields_vtk("waveguide_field.vtu", freq=10e9, active_port=0)
print("Saved: waveguide_field.vtu")
print("Open in ParaView to visualize E-field")
```

In ParaView:
1. Open `waveguide_field.vtu`
2. Apply "Surface" representation
3. Color by "E_field" magnitude
4. Add "Glyph" filter to show field vectors

## Step 8: Low-Level API (Optional)

For more control, use the core API directly:

```python
import pyedgefem as em

# Load existing mesh
mesh = em.load_gmsh("waveguide_mesh.msh")

# Set up boundary conditions
bc = em.build_edge_pec(mesh, 1)  # Tag 1 = PEC walls

# Configure port dimensions
port_dim = em.RectWaveguidePort()
port_dim.a = a
port_dim.b = b

# Solve TE10 mode at frequency
freq = 10e9
mode = em.solve_te10_mode(port_dim, freq)

# Build ports
surface1 = em.extract_surface_mesh(mesh, 2)  # Tag 2 = Port 1
em.populate_te10_field(surface1, port_dim, mode)
port1 = em.build_wave_port(mesh, surface1, mode)

surface2 = em.extract_surface_mesh(mesh, 3)  # Tag 3 = Port 2
em.populate_te10_field(surface2, port_dim, mode)
port2 = em.build_wave_port(mesh, surface2, mode)

# Maxwell parameters
params = em.MaxwellParams()
params.omega = 2 * np.pi * freq
params.port_abc_scale = 0.5  # Optimal for accuracy

# Compute S-parameters
S = em.calculate_sparams_eigenmode(mesh, params, bc, [port1, port2])
print(f"S21 = {S[1,0]:.6f}")
```

## Troubleshooting

### |S11| > 0 dB (non-physical)

- **Cause**: Mesh too coarse or port weights not matching field distribution
- **Solution**: Increase mesh density, use `use_eigenmode=True`

### |S21| << 1.0 (large insertion loss)

- **Cause**: Frequency below cutoff, or numerical issues
- **Solution**: Check `freq > wg.cutoff_frequency`, verify mesh quality

### Phase error > 10°

- **Cause**: Mesh too coarse for accurate phase
- **Solution**: Increase mesh density to 15-20 elements/wavelength

### Gmsh error

- **Cause**: Gmsh not installed or not in PATH
- **Solution**: Install Gmsh: `brew install gmsh` (macOS) or `apt install gmsh` (Linux)

## Complete Script

```python
"""
Complete WR-90 Waveguide Analysis Script
"""
from edgefem.designs import RectWaveguideDesign
import edgefem.plots as vp
import numpy as np
import matplotlib.pyplot as plt

# Design parameters
a = 22.86e-3  # 22.86 mm
b = 10.16e-3  # 10.16 mm
L = 50e-3     # 50 mm

# Create and analyze
wg = RectWaveguideDesign(a=a, b=b, length=L)
print(f"WR-90 Waveguide: {a*1000:.2f} x {b*1000:.2f} x {L*1000:.2f} mm")
print(f"TE10 cutoff: {wg.cutoff_frequency/1e9:.3f} GHz")

# Generate mesh
wg.generate_mesh(density=12)
print(f"Mesh: {len(wg.mesh.elements)} elements")

# Frequency sweep
freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=41, verbose=True)

# Plot
vp.plot_sparams_vs_freq(
    freqs, S_list,
    params=['S11', 'S21'],
    title=f"WR-90 Waveguide ({L*1000:.0f}mm)",
    save="wr90_complete.png"
)

# Export
wg.export_touchstone("wr90_complete", freqs, S_list)
wg.export_fields_vtk("wr90_field.vtu", freq=10e9)

print("\nAnalysis complete!")
print("Files saved: wr90_complete.s2p, wr90_complete.png, wr90_field.vtu")
```

## Next Steps

- [Patch Antenna Tutorial](patch_antenna.md) - Design antennas
- [Unit Cell Tutorial](unit_cell.md) - Periodic structures
- [API Reference](../api/core.md) - Full function documentation
