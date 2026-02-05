# Tutorial: Periodic Unit Cell Analysis

This tutorial covers analyzing periodic structures (metasurfaces, frequency selective surfaces, phased array unit cells) with VectorEM. You'll learn to:

- Create unit cell geometry with periodic boundaries
- Compute reflection and transmission coefficients
- Analyze oblique incidence behavior
- Extract surface impedance and effective parameters
- Study scan angle effects for phased arrays

## Background

Periodic structures are fundamental building blocks for:

- **Metasurfaces**: Control wave phase, polarization, absorption
- **Frequency Selective Surfaces (FSS)**: Bandpass/bandstop filters
- **Phased Arrays**: Unit cell characterizes element coupling
- **Metamaterials**: Engineered effective ε and μ

VectorEM uses Floquet boundary conditions to simulate an infinite periodic array from a single unit cell.

**Key concepts:**

- Floquet theorem: Fields repeat with phase shift across unit cell
- Phase shift: $\phi = k_0 d \sin\theta$ for incidence angle θ
- Normal incidence: θ = 0, no phase shift between cells

## Step 1: Create Unit Cell

Design a simple square patch unit cell for 10 GHz operation:

```python
from vectorem.designs import UnitCellDesign
import numpy as np

# Unit cell parameters
period = 5e-3        # 5mm period (λ/6 at 10 GHz)
h_sub = 0.5e-3       # 0.5mm substrate
eps_r = 3.5          # Low-loss dielectric

# Create unit cell with ground plane (reflective)
cell = UnitCellDesign(
    period_x=period,
    period_y=period,
    substrate_height=h_sub,
    substrate_eps_r=eps_r,
    has_ground_plane=True,  # Reflective metasurface
)

print(f"Unit cell: {period*1000:.1f} x {period*1000:.1f} mm")
print(f"Substrate: h={h_sub*1000:.2f}mm, εr={eps_r}")
```

## Step 2: Add Patch Element

Add a square patch to the unit cell:

```python
# Square patch (80% fill factor)
patch_size = 0.8 * period  # 4mm

cell.add_patch(
    width=patch_size,
    length=patch_size,
    center_x=0,    # Centered
    center_y=0,
    rotation=0,    # No rotation
    layer="top"    # On substrate top surface
)

print(f"Patch size: {patch_size*1000:.1f} x {patch_size*1000:.1f} mm")
```

## Step 3: Generate Mesh

```python
# Generate mesh with periodic boundaries
cell.generate_mesh(
    density=20,              # Elements per wavelength
    design_freq=10e9,        # Design frequency for mesh sizing
    output_path="unit_cell"  # Save mesh for inspection
)

print(f"Mesh generated: {len(cell.mesh.elements)} elements")
```

## Step 4: Normal Incidence Analysis

Compute reflection coefficient at normal incidence:

```python
# Normal incidence (θ=0, φ=0)
freq = 10e9
R, T = cell.reflection_transmission(freq, theta=0, phi=0, pol='TE')

# For grounded structure, T ≈ 0
print(f"\nNormal incidence at {freq/1e9:.1f} GHz:")
print(f"  Reflection: R = {R:.4f}")
print(f"  |R| = {abs(R):.4f} ({20*np.log10(abs(R)+1e-10):.1f} dB)")
print(f"  Phase(R) = {np.angle(R)*180/np.pi:.1f}°")
```

## Step 5: Frequency Sweep

Study frequency-dependent behavior:

```python
# Frequency sweep
freqs, R_array, T_array = cell.frequency_sweep_reflection(
    f_start=5e9,
    f_stop=15e9,
    n_points=51,
    theta=0,
    phi=0,
    pol='TE',
    verbose=True
)

# Find resonance (minimum reflection or phase crossing)
R_mag = np.abs(R_array)
idx_min = np.argmin(R_mag)
f_resonance = freqs[idx_min]

print(f"\nResonance: {f_resonance/1e9:.2f} GHz")
print(f"  |R|_min = {R_mag[idx_min]:.4f}")
```

Plot the results:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Magnitude
axes[0].plot(freqs/1e9, 20*np.log10(R_mag), 'b-', linewidth=2)
axes[0].set_xlabel('Frequency (GHz)')
axes[0].set_ylabel('|R| (dB)')
axes[0].set_title('Unit Cell Reflection Magnitude')
axes[0].grid(True)
axes[0].axvline(x=f_resonance/1e9, color='r', linestyle='--', label='Resonance')
axes[0].legend()

# Phase
R_phase = np.angle(R_array) * 180 / np.pi
axes[1].plot(freqs/1e9, R_phase, 'b-', linewidth=2)
axes[1].set_xlabel('Frequency (GHz)')
axes[1].set_ylabel('Phase(R) (deg)')
axes[1].set_title('Reflection Phase')
axes[1].grid(True)
axes[1].axhline(y=0, color='r', linestyle='--')
axes[1].axhline(y=-180, color='g', linestyle=':')

plt.tight_layout()
plt.savefig('unit_cell_reflection.png', dpi=150)
plt.show()
```

## Step 6: Surface Impedance

Extract surface impedance from reflection:

```python
# Surface impedance at resonance
Zs = cell.surface_impedance(f_resonance)

print(f"\nSurface impedance at {f_resonance/1e9:.2f} GHz:")
print(f"  Zs = {Zs.real:.1f} + j{Zs.imag:.1f} Ω")

# For reference, free-space impedance
Z0 = 376.73  # Ω
print(f"  Normalized: Zs/Z0 = {(Zs/Z0).real:.3f} + j{(Zs/Z0).imag:.3f}")
```

Surface impedance characterization:

| Zs Behavior | Metasurface Type |
|-------------|-----------------|
| Re(Zs) >> Z0 | High-impedance surface (HIS) |
| Re(Zs) << Z0 | Low-impedance (conductor-like) |
| Im(Zs) = 0 | Resonance condition |
| Im(Zs) > 0 | Inductive |
| Im(Zs) < 0 | Capacitive |

## Step 7: Oblique Incidence

Analyze behavior vs. incidence angle:

```python
# Scan angle sweep at resonant frequency
angles = np.arange(0, 61, 5)  # 0 to 60 degrees
R_vs_angle_TE = []
R_vs_angle_TM = []

for theta in angles:
    R_te, _ = cell.reflection_transmission(f_resonance, theta=theta, phi=0, pol='TE')
    R_tm, _ = cell.reflection_transmission(f_resonance, theta=theta, phi=0, pol='TM')
    R_vs_angle_TE.append(R_te)
    R_vs_angle_TM.append(R_tm)

R_vs_angle_TE = np.array(R_vs_angle_TE)
R_vs_angle_TM = np.array(R_vs_angle_TM)
```

Plot angular dependence:

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Magnitude vs angle
axes[0].plot(angles, np.abs(R_vs_angle_TE), 'b-o', label='TE')
axes[0].plot(angles, np.abs(R_vs_angle_TM), 'r--s', label='TM')
axes[0].set_xlabel('Incidence Angle (deg)')
axes[0].set_ylabel('|R|')
axes[0].set_title(f'Reflection Magnitude at {f_resonance/1e9:.1f} GHz')
axes[0].legend()
axes[0].grid(True)

# Phase vs angle
axes[1].plot(angles, np.angle(R_vs_angle_TE)*180/np.pi, 'b-o', label='TE')
axes[1].plot(angles, np.angle(R_vs_angle_TM)*180/np.pi, 'r--s', label='TM')
axes[1].set_xlabel('Incidence Angle (deg)')
axes[1].set_ylabel('Phase(R) (deg)')
axes[1].set_title('Reflection Phase')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig('unit_cell_oblique.png', dpi=150)
plt.show()
```

## Step 8: FSS Configuration (Transmission)

For frequency selective surfaces that allow transmission, configure without ground plane:

```python
# FSS unit cell (no ground plane)
fss = UnitCellDesign(
    period_x=5e-3,
    period_y=5e-3,
    substrate_height=0.5e-3,
    substrate_eps_r=3.5,
    has_ground_plane=False,      # No ground
    air_height_below=15e-3,      # Air region below
)

# Add patch
fss.add_patch(width=4e-3, length=4e-3)

# Generate mesh
fss.generate_mesh(density=20)

# Compute R and T
R, T = fss.reflection_transmission(10e9, theta=0, phi=0)

print(f"\nFSS at 10 GHz:")
print(f"  |R| = {abs(R):.4f} ({20*np.log10(abs(R)+1e-10):.1f} dB)")
print(f"  |T| = {abs(T):.4f} ({20*np.log10(abs(T)+1e-10):.1f} dB)")

# Energy conservation check
power = abs(R)**2 + abs(T)**2
print(f"  |R|² + |T|² = {power:.4f} (should be ≤1)")
```

## Step 9: Effective Parameter Extraction

For FSS structures, extract effective ε and μ using NRW method:

```python
if fss.is_fss:
    eps_eff, mu_eff = fss.effective_parameters(10e9)

    print(f"\nEffective parameters at 10 GHz:")
    print(f"  εeff = {eps_eff.real:.3f} + j{eps_eff.imag:.3f}")
    print(f"  μeff = {mu_eff.real:.3f} + j{mu_eff.imag:.3f}")

    # Refractive index
    n_eff = np.sqrt(eps_eff * mu_eff)
    print(f"  n_eff = {n_eff.real:.3f} + j{n_eff.imag:.3f}")
```

## Step 10: Export Results

### Touchstone Format

```python
# Export reflection data
cell.export_touchstone("unit_cell_reflection", freqs, R_array)
print("Saved: unit_cell_reflection.s1p")

# For FSS with transmission
fss.export_touchstone("fss_sparams", freqs, R_array, T_array)
print("Saved: fss_sparams.s2p")
```

## Advanced: Complex Unit Cell

Create a more complex unit cell with multiple elements:

```python
# Jerusalem cross unit cell
jc = UnitCellDesign(
    period_x=10e-3,
    period_y=10e-3,
    substrate_height=1.0e-3,
    substrate_eps_r=2.2,
)

# Central cross (horizontal arm)
jc.add_patch(width=8e-3, length=2e-3, center_x=0, center_y=0)

# Vertical arm
jc.add_patch(width=2e-3, length=8e-3, center_x=0, center_y=0)

# End caps (makes it a Jerusalem cross)
jc.add_patch(width=2e-3, length=4e-3, center_x=-3e-3, center_y=0)
jc.add_patch(width=2e-3, length=4e-3, center_x=3e-3, center_y=0)

jc.generate_mesh(density=20)
print(f"Jerusalem cross: {len(jc._geometry_elements)} elements")
```

## Phased Array Unit Cell Analysis

For phased array applications, analyze active impedance vs. scan:

```python
from vectorem.designs import UnitCellDesign
import numpy as np

# Phased array unit cell
array_cell = UnitCellDesign(
    period_x=15e-3,      # ~λ/2 at 10 GHz
    period_y=15e-3,
    substrate_height=1.5e-3,
    substrate_eps_r=3.0,
)

# Patch element
array_cell.add_patch(width=10e-3, length=10e-3)
array_cell.generate_mesh(density=15)

# Active impedance vs. E-plane scan
scan_angles = np.arange(0, 71, 5)
Z_active = []

for theta in scan_angles:
    Zs = array_cell.surface_impedance(10e9, theta=theta, phi=0, pol='TE')
    Z_active.append(Zs)

Z_active = np.array(Z_active)

# Plot active impedance
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(scan_angles, Z_active.real, 'b-o', label='Re(Z_active)')
ax.plot(scan_angles, Z_active.imag, 'r--s', label='Im(Z_active)')
ax.set_xlabel('Scan Angle (deg)')
ax.set_ylabel('Impedance (Ω)')
ax.set_title('Active Element Impedance vs. Scan Angle')
ax.legend()
ax.grid(True)
plt.savefig('active_impedance_scan.png', dpi=150)
plt.show()

# Find scan blindness (impedance anomaly)
Z_mag = np.abs(Z_active)
if np.max(Z_mag) > 5 * np.min(Z_mag):
    idx_blind = np.argmax(Z_mag)
    print(f"Potential scan blindness at θ = {scan_angles[idx_blind]}°")
```

## Design Guidelines

### Unit Cell Size

| Application | Period | Notes |
|-------------|--------|-------|
| Metasurface | < λ/5 | Subwavelength for effective medium |
| FSS | λ/3 - λ/2 | Resonant behavior |
| Phased array | ~λ/2 | Grating lobe avoidance |

### Element Geometry

| Type | Behavior | Application |
|------|----------|-------------|
| Patch | Resonant, capacitive below | Absorber, HIS |
| Slot | Resonant, inductive below | FSS, AMC |
| Cross | Dual-polarization | Polarizer |
| Ring | Multiple resonances | Multi-band |

### Substrate Effects

- Higher εr → narrower bandwidth, smaller cell
- Thicker substrate → stronger resonance, bandwidth
- Loss tangent → absorption, Q reduction

## Troubleshooting

### Non-physical results (|R| > 1)

- Increase mesh density
- Check periodic boundary matching
- Verify element fits within unit cell

### No resonance visible

- Element may be too small/large for frequency
- Check element is properly meshed
- Verify substrate properties

### TE/TM results identical

- At normal incidence, TE and TM should match for symmetric cells
- Check polarization definition for oblique incidence

## Complete Script

```python
"""
Complete Unit Cell Analysis Script
"""
from vectorem.designs import UnitCellDesign
import numpy as np
import matplotlib.pyplot as plt

# Create unit cell
cell = UnitCellDesign(
    period_x=5e-3,
    period_y=5e-3,
    substrate_height=0.5e-3,
    substrate_eps_r=3.5,
)

cell.add_patch(width=4e-3, length=4e-3)
cell.generate_mesh(density=20)

print(f"Unit cell: {cell.period_x*1000:.1f} x {cell.period_y*1000:.1f} mm")

# Frequency sweep
freqs, R_arr, _ = cell.frequency_sweep_reflection(
    5e9, 15e9, n_points=51, verbose=True
)

# Find resonance
idx_res = np.argmin(np.abs(R_arr))
f_res = freqs[idx_res]
print(f"\nResonance: {f_res/1e9:.2f} GHz")

# Surface impedance
Zs = cell.surface_impedance(f_res)
print(f"Surface impedance: {Zs.real:.1f}+j{Zs.imag:.1f} Ω")

# Oblique incidence
angles = np.arange(0, 61, 10)
R_vs_theta = [cell.reflection_transmission(f_res, theta=a)[0] for a in angles]

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Magnitude vs freq
axes[0,0].plot(freqs/1e9, 20*np.log10(np.abs(R_arr)+1e-10))
axes[0,0].set_xlabel('Frequency (GHz)')
axes[0,0].set_ylabel('|R| (dB)')
axes[0,0].set_title('Reflection vs Frequency')
axes[0,0].grid(True)

# Phase vs freq
axes[0,1].plot(freqs/1e9, np.angle(R_arr)*180/np.pi)
axes[0,1].set_xlabel('Frequency (GHz)')
axes[0,1].set_ylabel('Phase(R) (deg)')
axes[0,1].set_title('Phase vs Frequency')
axes[0,1].grid(True)

# |R| vs angle
axes[1,0].plot(angles, np.abs(R_vs_theta), 'o-')
axes[1,0].set_xlabel('Incidence Angle (deg)')
axes[1,0].set_ylabel('|R|')
axes[1,0].set_title(f'Reflection vs Angle at {f_res/1e9:.1f} GHz')
axes[1,0].grid(True)

# Phase vs angle
axes[1,1].plot(angles, np.angle(R_vs_theta)*180/np.pi, 'o-')
axes[1,1].set_xlabel('Incidence Angle (deg)')
axes[1,1].set_ylabel('Phase(R) (deg)')
axes[1,1].set_title('Phase vs Angle')
axes[1,1].grid(True)

plt.tight_layout()
plt.savefig('unit_cell_complete.png', dpi=150)
plt.show()

# Export
cell.export_touchstone("unit_cell", freqs, R_arr)
print("\nSaved: unit_cell.s1p, unit_cell_complete.png")
```

## Next Steps

- [Waveguide Tutorial](waveguide.md) - S-parameter analysis
- [Patch Antenna Tutorial](patch_antenna.md) - Antenna design
- [API Reference](../api/designs.md) - Full class documentation
