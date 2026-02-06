# Tutorial: Microstrip Patch Antenna Design

This tutorial covers designing and simulating microstrip patch antennas with EdgeFEM. You'll learn to:

- Design a patch antenna for a target frequency
- Configure feed mechanisms (probe and edge feed)
- Compute input impedance and return loss
- Generate radiation patterns
- Analyze bandwidth and directivity

## Background

Microstrip patch antennas are popular for their low profile, ease of fabrication, and integration with PCBs. A rectangular patch antenna resonates when its length is approximately half a wavelength in the dielectric medium.

**Design equations:**

- Effective dielectric constant: $\varepsilon_{eff} = \frac{\varepsilon_r + 1}{2} + \frac{\varepsilon_r - 1}{2}\left(1 + 12\frac{h}{W}\right)^{-1/2}$
- Fringing extension: $\Delta L = 0.412h \frac{(\varepsilon_{eff} + 0.3)(W/h + 0.264)}{(\varepsilon_{eff} - 0.258)(W/h + 0.8)}$
- Resonant frequency: $f_r = \frac{c}{2(L + 2\Delta L)\sqrt{\varepsilon_{eff}}}$

## Step 1: Design for Target Frequency

Let's design a 2.4 GHz patch antenna on FR-4 substrate:

```python
from edgefem.designs import PatchAntennaDesign
import numpy as np

# Target frequency
f_target = 2.4e9  # 2.4 GHz (WiFi band)

# Substrate parameters (FR-4)
eps_r = 4.4       # Relative permittivity
tan_d = 0.02      # Loss tangent
h = 1.6e-3        # Thickness: 1.6 mm (standard)

# Initial patch dimensions (will be refined)
# Rule of thumb: L ≈ c / (2 * f * sqrt(eps_eff))
c0 = 299792458.0
eps_eff_approx = (eps_r + 1) / 2
L_approx = c0 / (2 * f_target * np.sqrt(eps_eff_approx))
W_approx = 1.3 * L_approx  # Width slightly larger for better efficiency

print(f"Initial estimate: L = {L_approx*1000:.1f} mm, W = {W_approx*1000:.1f} mm")
```

Create the patch antenna:

```python
patch = PatchAntennaDesign(
    patch_length=29e-3,      # 29 mm (resonant dimension)
    patch_width=38e-3,       # 38 mm
    substrate_height=h,
    substrate_eps_r=eps_r,
    substrate_tan_d=tan_d,
)

# EdgeFEM estimates resonant frequency
print(f"Estimated resonant frequency: {patch.estimated_resonant_frequency/1e9:.3f} GHz")
```

## Step 2: Configure Feed

### Option A: Coaxial Probe Feed

The probe feed connects the ground plane to the patch through the substrate:

```python
# Probe at optimal position for ~50 ohm match
# Position is relative to patch center
patch.set_probe_feed(
    x_offset=0,          # Centered in width
    y_offset=-5e-3,      # 5mm from center toward edge
    probe_radius=0.5e-3  # 0.5mm radius probe
)
```

Feed position affects input impedance:
- At edge: highest impedance (~200-400 Ω)
- At center: zero impedance (null point)
- Intermediate: use for 50 Ω match

### Option B: Edge Feed with Inset

For microstrip feed line integration:

```python
# Alternative: inset edge feed
patch_edge = PatchAntennaDesign(
    patch_length=29e-3,
    patch_width=38e-3,
    substrate_height=h,
    substrate_eps_r=eps_r,
)

patch_edge.set_edge_feed(
    inset_length=8e-3,   # 8mm inset for impedance matching
    line_width=3e-3,     # 3mm feed line (≈50Ω on FR-4)
    edge='y_neg'         # Feed from -Y edge
)
```

## Step 3: Generate Mesh

```python
# Generate mesh (15 elements per wavelength recommended for antennas)
patch.generate_mesh(density=15)

print(f"Mesh generated: {len(patch.mesh.elements)} elements")
```

Save mesh for inspection:

```python
patch.generate_mesh(density=15, output_path="patch_antenna")
# Creates: patch_antenna.geo, patch_antenna.msh
```

## Step 4: Input Impedance Analysis

```python
# Compute input impedance at design frequency
Z_in = patch.input_impedance(2.4e9)

print(f"\nInput impedance at 2.4 GHz:")
print(f"  Z_in = {Z_in.real:.1f} + j{Z_in.imag:.1f} Ω")

# Reflection coefficient (assuming 50Ω system)
Z0 = 50
Gamma = (Z_in - Z0) / (Z_in + Z0)
print(f"\nReflection coefficient:")
print(f"  |Γ| = {abs(Gamma):.4f}")
print(f"  Return loss = {20*np.log10(abs(Gamma)+1e-10):.1f} dB")

# VSWR
vswr = (1 + abs(Gamma)) / (1 - abs(Gamma))
print(f"  VSWR = {vswr:.2f}")
```

## Step 5: Return Loss Sweep

Find the resonant frequency and bandwidth:

```python
# Frequency sweep around design frequency
freqs, Z_in_array, S11_array = patch.frequency_sweep(
    f_start=2.0e9,
    f_stop=2.8e9,
    n_points=81,
    z0=50.0,
    verbose=True
)

# Find minimum return loss (resonant frequency)
rl_db = 20 * np.log10(np.abs(S11_array) + 1e-10)
idx_min = np.argmin(rl_db)
f_resonant = freqs[idx_min]

print(f"\nResonant frequency: {f_resonant/1e9:.3f} GHz")
print(f"Minimum return loss: {rl_db[idx_min]:.1f} dB")
```

Plot the results:

```python
import matplotlib.pyplot as plt
import edgefem.plots as vp

# Return loss plot
vp.plot_return_loss(freqs, S11_array,
                    title="Patch Antenna Return Loss",
                    save="patch_return_loss.png")
```

Custom impedance plot:

```python
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Return loss
axes[0].plot(freqs/1e9, rl_db, 'b-', linewidth=2)
axes[0].axhline(y=-10, color='r', linestyle='--', label='-10 dB')
axes[0].set_xlabel('Frequency (GHz)')
axes[0].set_ylabel('Return Loss (dB)')
axes[0].set_title('Patch Antenna Return Loss')
axes[0].legend()
axes[0].grid(True)
axes[0].set_ylim(-30, 0)

# Input impedance
axes[1].plot(freqs/1e9, Z_in_array.real, 'b-', label='Re(Z)', linewidth=2)
axes[1].plot(freqs/1e9, Z_in_array.imag, 'r--', label='Im(Z)', linewidth=2)
axes[1].axhline(y=50, color='g', linestyle=':', label='50Ω')
axes[1].set_xlabel('Frequency (GHz)')
axes[1].set_ylabel('Impedance (Ω)')
axes[1].set_title('Input Impedance')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig('patch_impedance.png', dpi=150)
plt.show()
```

## Step 6: Bandwidth Analysis

```python
# Compute bandwidth for VSWR < 2 (return loss < -10 dB)
f_low, f_high, bw_percent = patch.bandwidth(
    vswr_threshold=2.0,
    z0=50.0
)

print(f"\nBandwidth (VSWR < 2.0):")
print(f"  Lower frequency: {f_low/1e9:.3f} GHz")
print(f"  Upper frequency: {f_high/1e9:.3f} GHz")
print(f"  Bandwidth: {bw_percent:.1f}%")
print(f"  Bandwidth: {(f_high - f_low)/1e6:.1f} MHz")
```

## Step 7: Radiation Pattern

Compute the 3D far-field pattern:

```python
# Compute radiation pattern at resonance
pattern = patch.radiation_pattern(f_resonant, n_theta=91, n_phi=73)

# Compute directivity
D = patch.directivity(f_resonant)
D_dBi = 10 * np.log10(D)

print(f"\nRadiation characteristics at {f_resonant/1e9:.3f} GHz:")
print(f"  Directivity: {D_dBi:.1f} dBi")
```

Plot pattern cuts:

```python
# E-plane and H-plane cuts
vp.plot_pattern_3d_cuts(
    pattern,
    title=f"Patch Antenna Pattern ({f_resonant/1e9:.2f} GHz)",
    min_db=-30,
    save="patch_pattern_cuts.png"
)
```

Interactive 3D pattern (requires plotly):

```python
try:
    vp.plot_pattern_3d(
        pattern,
        title="Patch Antenna 3D Pattern",
        interactive=True,
        save="patch_pattern_3d.html"
    )
    print("Interactive 3D pattern saved to patch_pattern_3d.html")
except ImportError:
    print("Install plotly for 3D interactive: pip install plotly")
```

## Step 8: Export Results

### Touchstone Format

```python
# Export S11 to Touchstone
patch.export_touchstone("patch_antenna", freqs, S11_array)
print("Saved: patch_antenna.s1p")
```

### Pattern Data

```python
# Export pattern to CSV
patch.export_pattern_csv("patch_pattern.csv", f_resonant)
print("Saved: patch_pattern.csv")
```

## Design Optimization

To adjust the design for better matching:

```python
# Parametric study: vary feed position
feed_positions = np.linspace(-10e-3, -3e-3, 8)
results = []

for y_pos in feed_positions:
    # Create new patch with different feed position
    p = PatchAntennaDesign(
        patch_length=29e-3,
        patch_width=38e-3,
        substrate_height=1.6e-3,
        substrate_eps_r=4.4,
    )
    p.set_probe_feed(x_offset=0, y_offset=y_pos)
    p.generate_mesh(density=12)

    Z_in = p.input_impedance(2.4e9)
    rl = p.return_loss(2.4e9)

    results.append({
        'y_pos': y_pos * 1000,  # mm
        'R': Z_in.real,
        'X': Z_in.imag,
        'RL': rl
    })
    print(f"y = {y_pos*1000:.1f} mm: Z = {Z_in.real:.0f}+j{Z_in.imag:.0f} Ω, RL = {rl:.1f} dB")

# Find best match
best = min(results, key=lambda r: abs(r['R'] - 50) + abs(r['X']))
print(f"\nBest match at y = {best['y_pos']:.1f} mm: Z = {best['R']:.0f}+j{best['X']:.0f} Ω")
```

## Design Guidelines

### Substrate Selection

| Substrate | εr | tan δ | Use Case |
|-----------|-----|-------|----------|
| FR-4 | 4.4 | 0.02 | Low cost, WiFi |
| Rogers RO4003C | 3.55 | 0.0027 | RF/microwave |
| Rogers RT/duroid 5880 | 2.2 | 0.0009 | High frequency |
| Alumina | 9.8 | 0.0001 | mmWave |

### Dimensions

| Parameter | Typical Range | Effect |
|-----------|---------------|--------|
| Length L | λ_eff/2 | Sets resonant frequency |
| Width W | 1.2-1.5 × L | Affects bandwidth, efficiency |
| Height h | 0.003-0.05 λ₀ | Higher = wider BW, more loss |

### Feed Position

| Feed Type | Impedance Range | Advantage |
|-----------|-----------------|-----------|
| Probe at edge | 200-400 Ω | Simple |
| Probe offset | 20-100 Ω | Direct 50Ω match |
| Inset edge | 20-100 Ω | Planar, PCB compatible |
| Aperture coupled | Wide range | Better BW, isolation |

## Troubleshooting

### Resonance too high/low

- **Too high**: Increase patch length
- **Too low**: Decrease patch length
- Rule: Δf/f ≈ -ΔL/L

### Can't achieve 50Ω match

- Adjust feed position (probe) or inset depth (edge feed)
- Consider different substrate thickness
- Use matching network

### Narrow bandwidth

- Increase substrate thickness
- Use lower εr substrate
- Consider stacked patch design

### Low directivity

- Increase patch width
- Check ground plane size (should be > λ/2 beyond patch)

## Complete Script

```python
"""
Complete Patch Antenna Design Script
"""
from edgefem.designs import PatchAntennaDesign
import edgefem.plots as vp
import numpy as np
import matplotlib.pyplot as plt

# Design parameters
f_design = 2.4e9  # Target frequency
eps_r = 4.4       # FR-4
h = 1.6e-3        # 1.6mm substrate

# Create antenna
patch = PatchAntennaDesign(
    patch_length=29e-3,
    patch_width=38e-3,
    substrate_height=h,
    substrate_eps_r=eps_r,
    substrate_tan_d=0.02,
)

print(f"Estimated resonance: {patch.estimated_resonant_frequency/1e9:.3f} GHz")

# Configure probe feed
patch.set_probe_feed(x_offset=0, y_offset=-5e-3)

# Generate mesh
patch.generate_mesh(density=15)
print(f"Mesh: {len(patch.mesh.elements)} elements")

# Frequency sweep
freqs, Z_in, S11 = patch.frequency_sweep(2.0e9, 2.8e9, n_points=81, verbose=True)

# Find resonance
rl_db = 20 * np.log10(np.abs(S11) + 1e-10)
f_res = freqs[np.argmin(rl_db)]
print(f"\nResonant frequency: {f_res/1e9:.3f} GHz")

# Bandwidth
f_lo, f_hi, bw = patch.bandwidth(vswr_threshold=2.0)
print(f"Bandwidth (VSWR<2): {bw:.1f}% ({(f_hi-f_lo)/1e6:.0f} MHz)")

# Pattern and directivity
pattern = patch.radiation_pattern(f_res)
D = patch.directivity(f_res)
print(f"Directivity: {10*np.log10(D):.1f} dBi")

# Plots
vp.plot_return_loss(freqs, S11, title="2.4 GHz Patch", save="patch_rl.png")
vp.plot_pattern_3d_cuts(pattern, title="Pattern Cuts", save="patch_cuts.png")

# Export
patch.export_touchstone("patch_2p4ghz", freqs, S11)

print("\nDesign complete!")
```

## Next Steps

- [Unit Cell Tutorial](unit_cell.md) - Periodic structures
- [Waveguide Tutorial](waveguide.md) - S-parameter analysis
- [Designs API](../api/designs.md) - Full class reference
