# Plotting Module API Reference

The `vectorem.plots` module provides publication-quality plotting functions for common EM visualization tasks.

```python
import vectorem.plots as vp
```

## Dependencies

- `matplotlib` (required)
- `numpy` (required)
- `plotly` (optional, for interactive 3D plots)

## Pattern Plots

### plot_pattern_2d

```python
vp.plot_pattern_2d(
    pattern,
    title: str = "Radiation Pattern",
    ylabel: str = "Gain (dB)",
    xlabel: str = "Angle (deg)",
    min_db: float = -40,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    **kwargs
)
```

Plot 2D radiation pattern in rectangular coordinates.

**Parameters:**

- `pattern` - Pattern object with `pattern_dB()` method or 1D array
- `title` - Plot title
- `ylabel` - Y-axis label
- `xlabel` - X-axis label
- `min_db` - Minimum dB value to display
- `figsize` - Figure size in inches
- `save` - Path to save figure (None to skip)
- `show` - Whether to display the figure
- `ax` - Matplotlib axes to plot on (creates new figure if None)

**Example:**
```python
import vectorem.plots as vp
import numpy as np

# Plot E-plane cut
theta = np.linspace(-90, 90, 181)
pattern_db = 10 * np.cos(np.deg2rad(theta))**2  # Example pattern
vp.plot_pattern_2d(pattern_db, title="E-plane", save="eplane.png")
```

---

### plot_pattern_polar

```python
vp.plot_pattern_polar(
    pattern,
    title: str = "Radiation Pattern",
    min_db: float = -40,
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional[plt.Axes] = None,
    **kwargs
)
```

Plot 2D radiation pattern in polar coordinates.

**Example:**
```python
# Polar plot of pattern
vp.plot_pattern_polar(pattern, title="H-plane", min_db=-30, save="hplane_polar.png")
```

---

### plot_pattern_3d

```python
vp.plot_pattern_3d(
    pattern: FFPattern3D,
    title: str = "3D Radiation Pattern",
    colorscale: str = "Viridis",
    min_db: float = -40,
    interactive: bool = True,
    save: Optional[str] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (10, 8)
)
```

Plot 3D radiation pattern as a sphere.

**Parameters:**

- `pattern` - `FFPattern3D` object from `stratton_chu_3d()`
- `title` - Plot title
- `colorscale` - Plotly colorscale name (for interactive mode)
- `min_db` - Minimum dB value for color scale
- `interactive` - Use Plotly for interactive 3D (requires plotly)
- `save` - Path to save figure
- `show` - Whether to display

**Example:**
```python
# 3D pattern from far-field computation
pattern = em.stratton_chu_3d(mesh, E_field, H_field, huygens_tag, freq)
vp.plot_pattern_3d(pattern, title="Patch Antenna", interactive=True)
```

---

### plot_pattern_3d_cuts

```python
vp.plot_pattern_3d_cuts(
    pattern: FFPattern3D,
    title: str = "Pattern Cuts",
    min_db: float = -40,
    figsize: Tuple[float, float] = (12, 5),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot E-plane and H-plane cuts of a 3D pattern side-by-side.

**Example:**
```python
# Plot principal plane cuts
vp.plot_pattern_3d_cuts(pattern, title="2.4 GHz Patch", save="pattern_cuts.png")
```

---

## S-Parameter Plots

### plot_sparams_vs_freq

```python
vp.plot_sparams_vs_freq(
    freqs: np.ndarray,
    S_matrices: List[np.ndarray],
    params: List[str] = ['S11', 'S21'],
    show_phase: bool = True,
    title: str = "S-Parameters",
    figsize: Tuple[float, float] = (10, 6),
    save: Optional[str] = None,
    show: bool = True,
    freq_unit: str = "GHz"
)
```

Plot S-parameter magnitude and phase vs. frequency.

**Parameters:**

- `freqs` - Frequency array in Hz
- `S_matrices` - List of S-matrices (one per frequency point)
- `params` - Which S-parameters to plot ('S11', 'S21', 'S12', 'S22')
- `show_phase` - Include phase subplot
- `title` - Plot title
- `freq_unit` - Frequency unit for axis label ('Hz', 'kHz', 'MHz', 'GHz')

**Example:**
```python
# Waveguide frequency sweep
wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
wg.generate_mesh(density=10)
freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=21)

vp.plot_sparams_vs_freq(
    freqs, S_list,
    params=['S11', 'S21'],
    title="WR-90 Waveguide",
    save="waveguide_sparams.png"
)
```

---

### plot_return_loss

```python
vp.plot_return_loss(
    freqs: np.ndarray,
    S11_array: np.ndarray,
    title: str = "Return Loss",
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot return loss (|S11| in dB) vs. frequency.

**Example:**
```python
# Plot antenna return loss
S11 = np.array([S[0,0] for S in S_list])
vp.plot_return_loss(freqs, S11, title="Patch Antenna", save="return_loss.png")
```

---

### plot_insertion_loss

```python
vp.plot_insertion_loss(
    freqs: np.ndarray,
    S21_array: np.ndarray,
    title: str = "Insertion Loss",
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot insertion loss (|S21| in dB) vs. frequency.

---

### plot_sparams_smith

```python
vp.plot_sparams_smith(
    S: np.ndarray,
    title: str = "Smith Chart",
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot S-parameters on a Smith chart.

!!! note
    This function is a stub for future implementation.

---

### plot_sparams_polar

```python
vp.plot_sparams_polar(
    S: np.ndarray,
    title: str = "S-Parameters (Polar)",
    figsize: Tuple[float, float] = (8, 8),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot S-parameters in polar format.

!!! note
    This function is a stub for future implementation.

---

## Coupling Plots

### plot_coupling_matrix

```python
vp.plot_coupling_matrix(
    S: np.ndarray,
    title: str = "Coupling Matrix",
    show_magnitude: bool = True,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    show: bool = True,
    annotate: bool = True
)
```

Plot S-parameter matrix as a heatmap.

**Parameters:**

- `S` - N×N S-parameter matrix
- `show_magnitude` - Show |S| if True, else show dB
- `cmap` - Matplotlib colormap
- `vmin`, `vmax` - Color scale limits
- `annotate` - Show values in cells

**Example:**
```python
# 4-element array coupling matrix
S = em.compute_coupling(mesh, params, bc, ports)
vp.plot_coupling_matrix(S, title="4-Element Array", save="coupling.png")
```

---

### plot_coupling_vs_distance

```python
vp.plot_coupling_vs_distance(
    coupling_data: np.ndarray,
    distances: np.ndarray,
    title: str = "Coupling vs. Distance",
    figsize: Tuple[float, float] = (8, 5),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot coupling coefficient vs. element separation.

**Example:**
```python
# Plot mutual coupling decay
distances = np.array([0.5, 1.0, 1.5, 2.0]) * wavelength
couplings = np.array([-15, -22, -28, -33])  # dB

vp.plot_coupling_vs_distance(couplings, distances,
                             title="Mutual Coupling Decay",
                             save="coupling_decay.png")
```

---

### plot_active_impedance_scan

```python
vp.plot_active_impedance_scan(
    scan_angles: np.ndarray,
    Z_active: np.ndarray,
    title: str = "Active Impedance vs. Scan",
    figsize: Tuple[float, float] = (10, 5),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot active impedance (real and imaginary) vs. scan angle.

**Example:**
```python
# Active impedance scan
angles = np.arange(0, 61, 5)
Z_active = em.active_impedance_scan(S, angles)

vp.plot_active_impedance_scan(angles, Z_active,
                              title="E-plane Scan",
                              save="active_z.png")
```

---

## Field Plots

### plot_field_slice

```python
vp.plot_field_slice(
    x: np.ndarray,
    y: np.ndarray,
    field: np.ndarray,
    title: str = "Field",
    xlabel: str = "x (mm)",
    ylabel: str = "y (mm)",
    cmap: str = "RdBu_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    colorbar_label: str = "Field (V/m)",
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    show: bool = True,
    ax: Optional[plt.Axes] = None
)
```

Plot a 2D field slice as a color map.

**Parameters:**

- `x`, `y` - Coordinate arrays (1D)
- `field` - 2D field array
- `cmap` - Colormap (default: diverging red-blue)
- `vmin`, `vmax` - Color scale limits

**Example:**
```python
# Plot E_z on xy-plane
vp.plot_field_slice(
    x_mm, y_mm, E_z_slice,
    title="E_z at z=0",
    cmap="RdBu_r",
    colorbar_label="E_z (V/m)",
    save="ez_slice.png"
)
```

---

### plot_field_magnitude

```python
vp.plot_field_magnitude(
    x: np.ndarray,
    y: np.ndarray,
    E_field: np.ndarray,
    title: str = "Electric Field Magnitude",
    figsize: Tuple[float, float] = (8, 6),
    save: Optional[str] = None,
    show: bool = True
)
```

Plot electric field magnitude |E| as a color map.

**Example:**
```python
# Plot |E| magnitude
E_mag = np.sqrt(np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2)
vp.plot_field_magnitude(x_mm, y_mm, E_mag,
                        title="E-field Magnitude",
                        save="e_magnitude.png")
```

---

## Complete Workflow Example

```python
import numpy as np
from vectorem.designs import RectWaveguideDesign, PatchAntennaDesign
import vectorem.plots as vp

# --- Waveguide S-parameters ---
wg = RectWaveguideDesign(a=22.86e-3, b=10.16e-3, length=50e-3)
wg.generate_mesh(density=10)
freqs, S_list = wg.frequency_sweep(8e9, 12e9, n_points=41)

# Multi-panel S-param plot
vp.plot_sparams_vs_freq(
    freqs, S_list,
    params=['S11', 'S21'],
    show_phase=True,
    title="WR-90 Waveguide S-Parameters",
    save="waveguide_analysis.png"
)

# --- Patch antenna pattern ---
patch = PatchAntennaDesign(
    patch_length=29e-3,
    patch_width=38e-3,
    substrate_height=1.6e-3,
    substrate_eps_r=4.4,
)
patch.set_probe_feed(y_offset=-5e-3)
patch.generate_mesh(density=15)

pattern = patch.radiation_pattern(2.4e9, n_theta=181, n_phi=361)

# Publication-quality pattern cuts
vp.plot_pattern_3d_cuts(
    pattern,
    title="2.4 GHz Microstrip Patch",
    min_db=-30,
    save="patch_pattern_cuts.png"
)

# Interactive 3D pattern (requires plotly)
try:
    vp.plot_pattern_3d(pattern, interactive=True, save="patch_3d.html")
except ImportError:
    print("Install plotly for interactive 3D: pip install plotly")
```
