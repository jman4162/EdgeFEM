# Design Classes API Reference

EdgeFEM provides high-level design classes that encapsulate common simulation workflows. These classes handle mesh generation, boundary conditions, and post-processing automatically.

## RectWaveguideDesign

High-level interface for rectangular waveguide S-parameter analysis.

```python
from edgefem.designs import RectWaveguideDesign
```

### Constructor

```python
RectWaveguideDesign(
    a: float,
    b: float,
    length: float,
    *,
    pec_tag: int = 1,
    port1_tag: int = 2,
    port2_tag: int = 3,
    volume_tag: int = 100
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `a` | float | Waveguide width (broadwall dimension) in meters |
| `b` | float | Waveguide height (narrowwall dimension) in meters |
| `length` | float | Waveguide section length in meters |
| `pec_tag` | int | Physical tag for PEC walls (default: 1) |
| `port1_tag` | int | Physical tag for port 1 surface (default: 2) |
| `port2_tag` | int | Physical tag for port 2 surface (default: 3) |
| `volume_tag` | int | Physical tag for waveguide volume (default: 100) |

**Standard Waveguide Dimensions:**

| Designation | Band | a (mm) | b (mm) | fc (GHz) |
|-------------|------|--------|--------|----------|
| WR-90 | X | 22.86 | 10.16 | 6.56 |
| WR-62 | Ku | 15.80 | 7.90 | 9.49 |
| WR-42 | Ka | 10.67 | 4.32 | 14.1 |
| WR-28 | Ka | 7.11 | 3.56 | 21.1 |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `mesh` | Mesh | FEM mesh (None until `generate_mesh` is called) |
| `ports` | List[WavePort] | Wave port definitions |
| `cutoff_frequency` | float | TE10 cutoff frequency in Hz |

### Methods

#### generate_mesh

```python
wg.generate_mesh(
    density: float = 10.0,
    output_path: Optional[str] = None,
    gmsh_options: Optional[List[str]] = None
)
```

Generate FEM mesh using Gmsh.

**Parameters:**

- `density` - Elements per wavelength at 1.5× cutoff frequency
- `output_path` - Path for mesh file (uses temp directory if None)
- `gmsh_options` - Additional Gmsh command-line options

#### load_mesh

```python
wg.load_mesh(mesh_path: str)
```

Load a pre-existing Gmsh v2 mesh file.

#### sparams_at_freq

```python
S = wg.sparams_at_freq(
    freq: float,
    *,
    use_eigenmode: bool = True,
    port_abc_scale: float = 0.5
) -> np.ndarray
```

Compute S-parameters at a single frequency.

**Parameters:**

- `freq` - Frequency in Hz (must be above cutoff)
- `use_eigenmode` - Use eigenmode-based extraction (recommended)
- `port_abc_scale` - ABC scaling factor (0.5 is optimal)

**Returns:** 2×2 complex S-parameter matrix

**Example:**
```python
S = wg.sparams_at_freq(10e9)
print(f"|S11| = {abs(S[0,0]):.4f}")
print(f"|S21| = {abs(S[1,0]):.4f}")
```

#### frequency_sweep

```python
freqs, S_list = wg.frequency_sweep(
    f_start: float,
    f_stop: float,
    n_points: int = 11,
    *,
    use_eigenmode: bool = True,
    port_abc_scale: float = 0.5,
    verbose: bool = False
) -> Tuple[np.ndarray, List[np.ndarray]]
```

Compute S-parameters over a frequency range.

**Returns:** Tuple of (frequencies array, list of S-matrices)

#### export_touchstone

```python
wg.export_touchstone(
    path: str,
    freqs: Optional[np.ndarray] = None,
    S_matrices: Optional[List[np.ndarray]] = None
)
```

Export S-parameters to Touchstone format (.s2p).

#### export_fields_vtk

```python
wg.export_fields_vtk(
    path: str,
    freq: Optional[float] = None,
    *,
    active_port: int = 0
)
```

Export E-field solution to VTK format for ParaView.

### Complete Example

```python
from edgefem.designs import RectWaveguideDesign
import edgefem.plots as vp
import numpy as np

# Create WR-90 waveguide
wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
print(f"Cutoff: {wg.cutoff_frequency/1e9:.2f} GHz")

# Generate mesh
wg.generate_mesh(density=10)

# Frequency sweep
freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=21, verbose=True)

# Plot results
vp.plot_sparams_vs_freq(freqs, S_list, params=['S11', 'S21'],
                        title="WR-90 Waveguide", save="wr90.png")

# Export
wg.export_touchstone("wr90_sweep", freqs, S_list)
wg.export_fields_vtk("wr90_field.vtu", freq=10e9)
```

---

## PatchAntennaDesign

High-level interface for microstrip patch antenna simulation.

```python
from edgefem.designs import PatchAntennaDesign
```

### Constructor

```python
PatchAntennaDesign(
    patch_length: float,
    patch_width: float,
    substrate_height: float,
    substrate_eps_r: float = 4.4,
    substrate_tan_d: float = 0.02,
    *,
    ground_plane_size: Optional[float] = None,
    air_box_height: Optional[float] = None
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `patch_length` | float | Patch length (resonant dimension) in meters |
| `patch_width` | float | Patch width in meters |
| `substrate_height` | float | Substrate thickness in meters |
| `substrate_eps_r` | float | Substrate relative permittivity (default: 4.4 for FR-4) |
| `substrate_tan_d` | float | Substrate loss tangent (default: 0.02) |
| `ground_plane_size` | float | Ground plane edge length (default: 3× patch) |
| `air_box_height` | float | Air box height above substrate |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `mesh` | Mesh | FEM mesh |
| `estimated_resonant_frequency` | float | Estimated resonance from transmission line model |

### Feed Configuration

#### set_probe_feed

```python
patch.set_probe_feed(
    x_offset: float = 0.0,
    y_offset: Optional[float] = None,
    probe_radius: float = 0.5e-3
)
```

Configure coaxial probe feed. Position is relative to patch center.

**Parameters:**

- `x_offset` - X offset from patch center (meters)
- `y_offset` - Y offset from patch center (None = optimal ~L/3 from edge)
- `probe_radius` - Probe radius (meters)

#### set_edge_feed

```python
patch.set_edge_feed(
    inset_length: float = 0.0,
    line_width: float = 3e-3,
    edge: str = 'y_neg'
)
```

Configure microstrip edge feed with optional inset.

**Parameters:**

- `inset_length` - Inset distance into patch (meters)
- `line_width` - Feed line width (meters)
- `edge` - Feed edge ('y_neg', 'y_pos', 'x_neg', 'x_pos')

### Simulation Methods

#### input_impedance

```python
Z_in = patch.input_impedance(freq: float) -> complex
```

Compute input impedance at given frequency (ohms).

#### return_loss

```python
rl_db = patch.return_loss(freq: float, z0: float = 50.0) -> float
```

Compute return loss |S11| in dB.

#### frequency_sweep

```python
freqs, Z_in, S11 = patch.frequency_sweep(
    f_start: float,
    f_stop: float,
    n_points: int = 51,
    z0: float = 50.0,
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

Sweep frequency and compute impedance and return loss.

#### radiation_pattern

```python
pattern = patch.radiation_pattern(
    freq: float,
    n_theta: int = 91,
    n_phi: int = 73
) -> FFPattern3D
```

Compute 3D far-field radiation pattern.

#### directivity

```python
D = patch.directivity(freq: float) -> float
```

Compute directivity (linear scale).

#### bandwidth

```python
f_low, f_high, bw_percent = patch.bandwidth(
    vswr_threshold: float = 2.0,
    z0: float = 50.0,
    search_range: Optional[Tuple[float, float]] = None
) -> Tuple[float, float, float]
```

Estimate bandwidth based on VSWR threshold.

### Complete Example

```python
from edgefem.designs import PatchAntennaDesign
import edgefem.plots as vp
import numpy as np

# Design 2.4 GHz patch on FR-4
patch = PatchAntennaDesign(
    patch_length=29e-3,
    patch_width=38e-3,
    substrate_height=1.6e-3,
    substrate_eps_r=4.4,
    substrate_tan_d=0.02,
)

print(f"Estimated resonance: {patch.estimated_resonant_frequency/1e9:.3f} GHz")

# Configure feed
patch.set_probe_feed(x_offset=0, y_offset=-5e-3)

# Generate mesh
patch.generate_mesh(density=15)

# Impedance analysis
Z_in = patch.input_impedance(2.4e9)
print(f"Z_in @ 2.4 GHz: {Z_in.real:.1f} + j{Z_in.imag:.1f} ohms")

rl = patch.return_loss(2.4e9)
print(f"Return loss: {rl:.1f} dB")

# Bandwidth
f_lo, f_hi, bw = patch.bandwidth(vswr_threshold=2.0)
print(f"Bandwidth (VSWR<2): {bw:.1f}%")

# Radiation pattern
pattern = patch.radiation_pattern(2.4e9)
D = patch.directivity(2.4e9)
print(f"Directivity: {10*np.log10(D):.1f} dBi")

# Plot pattern
vp.plot_pattern_3d_cuts(pattern, title="2.4 GHz Patch", save="patch_pattern.png")
```

---

## UnitCellDesign

High-level interface for periodic metasurface/metamaterial unit cell analysis.

```python
from edgefem.designs import UnitCellDesign
```

### Constructor

```python
UnitCellDesign(
    period_x: float,
    period_y: float,
    substrate_height: float = 1e-3,
    substrate_eps_r: float = 3.5,
    substrate_tan_d: float = 0.0,
    *,
    air_height_above: Optional[float] = None,
    air_height_below: Optional[float] = None,
    has_ground_plane: bool = True
)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `period_x` | float | Unit cell period in x-direction (meters) |
| `period_y` | float | Unit cell period in y-direction (meters) |
| `substrate_height` | float | Substrate thickness (meters) |
| `substrate_eps_r` | float | Substrate relative permittivity |
| `substrate_tan_d` | float | Substrate loss tangent |
| `air_height_above` | float | Air region above substrate |
| `air_height_below` | float | Air region below (for FSS configurations) |
| `has_ground_plane` | bool | Whether there's a PEC ground at z=0 |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `mesh` | Mesh | FEM mesh |
| `is_fss` | bool | True if FSS configuration (has transmission) |

### Geometry Methods

#### add_patch

```python
cell.add_patch(
    width: float,
    length: float,
    *,
    center_x: float = 0.0,
    center_y: float = 0.0,
    rotation: float = 0.0,
    layer: str = "top"
)
```

Add a rectangular PEC patch.

#### add_slot

```python
cell.add_slot(
    width: float,
    length: float,
    *,
    center_x: float = 0.0,
    center_y: float = 0.0,
    rotation: float = 0.0,
    layer: str = "ground"
)
```

Add a rectangular slot (aperture).

#### add_via

```python
cell.add_via(
    radius: float,
    *,
    center_x: float = 0.0,
    center_y: float = 0.0,
    top_layer: str = "top",
    bottom_layer: str = "ground"
)
```

Add a cylindrical via between layers.

### Simulation Methods

#### reflection_transmission

```python
R, T = cell.reflection_transmission(
    freq: float,
    *,
    theta: float = 0.0,
    phi: float = 0.0,
    pol: str = 'TE'
) -> Tuple[complex, complex]
```

Compute reflection and transmission coefficients.

**Parameters:**

- `freq` - Frequency in Hz
- `theta` - Incidence angle from normal (degrees)
- `phi` - Azimuth angle (degrees)
- `pol` - Polarization ('TE' or 'TM')

**Returns:** Tuple of (R, T) complex coefficients

#### frequency_sweep_reflection

```python
freqs, R_arr, T_arr = cell.frequency_sweep_reflection(
    f_start: float,
    f_stop: float,
    n_points: int = 51,
    *,
    theta: float = 0.0,
    phi: float = 0.0,
    pol: str = 'TE',
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

Sweep frequency and compute R/T coefficients.

#### surface_impedance

```python
Zs = cell.surface_impedance(freq: float, **kwargs) -> complex
```

Extract surface impedance from reflection coefficient (ohms).

Uses: Zs = Z0 × (1 + R) / (1 - R)

#### effective_parameters

```python
eps_eff, mu_eff = cell.effective_parameters(freq: float, **kwargs) -> Tuple[complex, complex]
```

Extract effective permittivity and permeability using NRW method.

!!! note
    Only valid for FSS configurations with transmission (`air_height_below` set, `has_ground_plane=False`).

### Complete Example

```python
from edgefem.designs import UnitCellDesign
import edgefem.plots as vp
import numpy as np

# Create unit cell for phased array
cell = UnitCellDesign(
    period_x=5e-3,
    period_y=5e-3,
    substrate_height=0.5e-3,
    substrate_eps_r=3.5,
)

# Add square patch
cell.add_patch(width=4e-3, length=4e-3)

# Generate mesh
cell.generate_mesh(density=20)

# Normal incidence reflection
R, T = cell.reflection_transmission(10e9, theta=0, phi=0)
print(f"|R| @ 10 GHz: {abs(R):.4f} ({20*np.log10(abs(R)+1e-10):.1f} dB)")

# Frequency sweep
freqs, R_arr, T_arr = cell.frequency_sweep_reflection(
    5e9, 15e9, n_points=51, verbose=True
)

# Plot reflection magnitude
import matplotlib.pyplot as plt
plt.figure()
plt.plot(freqs/1e9, 20*np.log10(np.abs(R_arr)+1e-10))
plt.xlabel('Frequency (GHz)')
plt.ylabel('|R| (dB)')
plt.title('Unit Cell Reflection')
plt.grid(True)
plt.savefig('unit_cell_reflection.png')

# Surface impedance
Zs = cell.surface_impedance(10e9)
print(f"Surface impedance: {Zs.real:.1f} + j{Zs.imag:.1f} ohms")

# Export to Touchstone
cell.export_touchstone("unit_cell", freqs, R_arr)
```

### Oblique Incidence Example

```python
# Scan angle sweep at 10 GHz
angles = np.arange(0, 61, 5)
R_vs_theta = []

for theta in angles:
    R, T = cell.reflection_transmission(10e9, theta=theta, phi=0, pol='TE')
    R_vs_theta.append(abs(R))

plt.figure()
plt.plot(angles, 20*np.log10(np.array(R_vs_theta)+1e-10))
plt.xlabel('Incidence Angle (deg)')
plt.ylabel('|R| (dB)')
plt.title('Reflection vs. Scan Angle')
plt.grid(True)
plt.savefig('scan_angle.png')
```
