# Core API Reference

This document covers the core functions and classes in the `pyedgefem` module.

## Mesh Operations

### load_gmsh

```python
mesh = em.load_gmsh(path: str) -> Mesh
```

Load a Gmsh v2 format mesh file.

**Parameters:**

- `path` - Path to `.msh` file (Gmsh v2 ASCII format)

**Returns:** `Mesh` object with nodes, elements, and edge connectivity

**Example:**
```python
mesh = em.load_gmsh("waveguide.msh")
print(f"Nodes: {len(mesh.nodes)}, Elements: {len(mesh.elements)}")
```

---

### validate_mesh

```python
report = em.validate_mesh(mesh: Mesh) -> MeshQualityReport
```

Check mesh quality and identify problematic elements.

**Parameters:**

- `mesh` - Mesh object to validate

**Returns:** `MeshQualityReport` with quality metrics

**Example:**
```python
report = em.validate_mesh(mesh)
if not report.is_valid():
    print(f"Warning: {report.num_degenerate_tets} degenerate elements")
```

---

### extract_surface_mesh

```python
surface = em.extract_surface_mesh(mesh: Mesh, tag: int) -> SurfaceMesh
```

Extract triangular surface mesh from a physical surface tag.

**Parameters:**

- `mesh` - Source mesh
- `tag` - Physical tag identifying the surface

**Returns:** `SurfaceMesh` containing surface triangles and edges

---

## Boundary Conditions

### build_edge_pec

```python
bc = em.build_edge_pec(mesh: Mesh, tag: int) -> BC
```

Build PEC (Perfect Electric Conductor) boundary condition for edges on a tagged surface.

**Parameters:**

- `mesh` - Source mesh
- `tag` - Physical tag identifying PEC surfaces

**Returns:** `BC` object with Dirichlet edge constraints

**Example:**
```python
# Tag 1 is PEC walls
bc = em.build_edge_pec(mesh, 1)
print(f"PEC edges: {len(bc.dirichlet_edges)}")
```

---

### BC class

```python
class BC:
    dirichlet_edges: Set[int]  # Edge indices with E_t = 0
```

Boundary condition container. Multiple BCs can be combined by merging their edge sets.

---

## Maxwell Assembly

### MaxwellParams

```python
class MaxwellParams:
    omega: float                    # Angular frequency (rad/s)
    eps_r_regions: Dict[int, complex]  # Relative permittivity per volume tag
    mu_r_regions: Dict[int, complex]   # Relative permeability per volume tag
    use_abc: bool                   # Enable ABC on outer boundaries
    use_pml: bool                   # Enable PML layers
    port_abc_scale: float           # ABC scaling factor (optimal: 0.5)
```

Maxwell equation parameters for assembly.

**Example:**
```python
params = em.MaxwellParams()
params.omega = 2 * np.pi * 10e9  # 10 GHz
params.eps_r_regions = {
    100: complex(1.0, 0),      # Air
    101: complex(4.4, -0.088), # FR-4 substrate (eps_r=4.4, tan_d=0.02)
}
params.port_abc_scale = 0.5  # Optimal for accurate S-parameters
```

---

### assemble_maxwell

```python
assembly = em.assemble_maxwell(
    mesh: Mesh,
    params: MaxwellParams,
    bc: BC,
    ports: List[WavePort],
    active_port: int = 0
) -> MaxwellAssembly
```

Assemble the Maxwell curl-curl FEM system matrix.

**Parameters:**

- `mesh` - FEM mesh
- `params` - Maxwell parameters
- `bc` - Boundary conditions
- `ports` - List of wave ports
- `active_port` - Index of port to excite (default: 0)

**Returns:** `MaxwellAssembly` with system matrix A and RHS vector b

---

## S-Parameter Calculation

### calculate_sparams_eigenmode

```python
S = em.calculate_sparams_eigenmode(
    mesh: Mesh,
    params: MaxwellParams,
    bc: BC,
    ports: List[WavePort]
) -> np.ndarray
```

Compute S-parameter matrix using eigenmode-based extraction. **Recommended** for highest accuracy.

**Parameters:**

- `mesh` - FEM mesh
- `params` - Maxwell parameters (set `port_abc_scale=0.5` for optimal accuracy)
- `bc` - Boundary conditions
- `ports` - List of wave ports

**Returns:** Complex S-matrix as `(N, N)` numpy array where N is number of ports

**Example:**
```python
params = em.MaxwellParams()
params.omega = 2 * np.pi * 10e9
params.port_abc_scale = 0.5  # Optimal

S = em.calculate_sparams_eigenmode(mesh, params, bc, [port1, port2])
print(f"S11 = {S[0,0]:.4f}, S21 = {S[1,0]:.4f}")
```

!!! note "Accuracy"
    This function achieves 99.4% transmission accuracy for matched waveguide sections.
    Always use `port_abc_scale=0.5` for best results.

---

### calculate_sparams

```python
S = em.calculate_sparams(
    mesh: Mesh,
    params: MaxwellParams,
    bc: BC,
    ports: List[WavePort]
) -> np.ndarray
```

Compute S-parameters using standard port extraction. Use `calculate_sparams_eigenmode` for higher accuracy.

---

### calculate_sparams_periodic

```python
S = em.calculate_sparams_periodic(
    mesh: Mesh,
    params: MaxwellParams,
    bc: BC,
    pbc: PeriodicBC,
    ports: List[WavePort]
) -> np.ndarray
```

Compute S-parameters for periodic structures with Floquet boundary conditions.

---

## Wave Ports

### RectWaveguidePort

```python
class RectWaveguidePort:
    a: float  # Waveguide width (m)
    b: float  # Waveguide height (m)
```

Rectangular waveguide port dimensions.

---

### solve_te10_mode

```python
mode = em.solve_te10_mode(port: RectWaveguidePort, freq: float) -> PortMode
```

Solve for TE10 mode parameters at given frequency.

**Parameters:**

- `port` - Waveguide dimensions
- `freq` - Frequency in Hz

**Returns:** `PortMode` with propagation constant, impedance, and field distribution

---

### build_wave_port

```python
port = em.build_wave_port(
    mesh: Mesh,
    surface: SurfaceMesh,
    mode: PortMode
) -> WavePort
```

Build a wave port from surface mesh and mode definition.

---

### populate_te10_field

```python
em.populate_te10_field(
    surface: SurfaceMesh,
    port: RectWaveguidePort,
    mode: PortMode
)
```

Populate TE10 field values on port surface edges.

---

## Linear Solver

### SolveOptions

```python
class SolveOptions:
    tolerance: float     # Convergence tolerance (default: 1e-8)
    max_iterations: int  # Maximum iterations (default: 10000)
    preconditioner: str  # "ILUT" or "diagonal"
```

---

### solve_linear

```python
result = em.solve_linear(
    A: SparseMatrix,
    b: Vector,
    options: Optional[SolveOptions] = None
) -> SolveResult
```

Solve the linear system Ax = b using BiCGSTAB with ILUT preconditioner.

**Returns:** `SolveResult` with solution vector `x`, convergence info

---

## Periodic Boundary Conditions

### build_periodic_pairs

```python
pbc = em.build_periodic_pairs(
    mesh: Mesh,
    tag_neg: int,
    tag_pos: int,
    period_vector: np.ndarray
) -> PeriodicBC
```

Build periodic boundary condition pairs between two surfaces.

**Parameters:**

- `mesh` - Source mesh
- `tag_neg` - Physical tag for -X (or -Y) surface
- `tag_pos` - Physical tag for +X (or +Y) surface
- `period_vector` - Translation vector [dx, dy, dz]

**Example:**
```python
pbc_x = em.build_periodic_pairs(mesh, 5, 6, np.array([5e-3, 0, 0]))
```

---

### set_floquet_phase

```python
em.set_floquet_phase(pbc: PeriodicBC, k_transverse: np.ndarray)
```

Set Floquet phase shift for oblique incidence.

**Parameters:**

- `pbc` - Periodic boundary condition
- `k_transverse` - Transverse wave vector [kx, ky, 0]

**Example:**
```python
# 30 degree incidence at 10 GHz
k0 = 2 * np.pi * 10e9 / 3e8
kx = k0 * np.sin(np.deg2rad(30))
em.set_floquet_phase(pbc_x, np.array([kx, 0, 0]))
```

---

## Far-Field Computation

### stratton_chu_3d

```python
pattern = em.stratton_chu_3d(
    mesh: Mesh,
    E_field: Vector,
    H_field: Vector,
    huygens_tag: int,
    freq: float,
    n_theta: int = 91,
    n_phi: int = 73
) -> FFPattern3D
```

Compute 3D far-field pattern using Stratton-Chu integration over a Huygens surface.

**Parameters:**

- `mesh` - FEM mesh
- `E_field` - Electric field solution
- `H_field` - Magnetic field solution (computed from E)
- `huygens_tag` - Physical tag of Huygens surface
- `freq` - Frequency in Hz
- `n_theta` - Number of theta samples (0 to 180 deg)
- `n_phi` - Number of phi samples (0 to 360 deg)

**Returns:** `FFPattern3D` with E_theta, E_phi components on spherical grid

---

### compute_directivity

```python
D = em.compute_directivity(pattern: FFPattern3D) -> float
```

Compute directivity from far-field pattern (linear scale).

---

### compute_max_gain

```python
G = em.compute_max_gain(pattern: FFPattern3D, efficiency: float = 1.0) -> float
```

Compute maximum gain (linear scale).

---

### compute_hpbw

```python
hpbw = em.compute_hpbw(pattern: FFPattern3D, plane: str = "E") -> float
```

Compute half-power beamwidth in degrees.

**Parameters:**

- `pattern` - Far-field pattern
- `plane` - "E" for E-plane or "H" for H-plane

---

## Export Functions

### write_touchstone_nport

```python
em.write_touchstone_nport(
    path: str,
    frequencies: List[float],
    S_matrices: List[np.ndarray]
)
```

Export S-parameters to Touchstone format.

**Parameters:**

- `path` - Output file path (.s1p, .s2p, etc.)
- `frequencies` - List of frequencies in Hz
- `S_matrices` - List of S-matrices (one per frequency)

---

### export_fields_vtk

```python
em.export_fields_vtk(
    mesh: Mesh,
    E_field: Vector,
    path: str,
    options: Optional[VTKExportOptions] = None
)
```

Export E-field solution to VTK format for ParaView visualization.

---

### write_pattern_3d_csv

```python
em.write_pattern_3d_csv(path: str, pattern: FFPattern3D)
```

Export 3D far-field pattern to CSV format.

---

### write_pattern_3d_vtk

```python
em.write_pattern_3d_vtk(path: str, pattern: FFPattern3D)
```

Export 3D far-field pattern to VTK format.

---

## Utility Functions

### check_passivity

```python
em.check_passivity(S: np.ndarray, tolerance: float = 0.01)
```

Verify that S-matrix satisfies passivity: |S11|² + |S21|² ≤ 1.

Raises `ValueError` if passivity is violated by more than tolerance.

---

### is_reciprocal

```python
reciprocal = em.is_reciprocal(S: np.ndarray, tolerance: float = 1e-3) -> bool
```

Check if S-matrix is reciprocal (S_ij = S_ji).

---

### to_db

```python
db_value = em.to_db(linear: float) -> float
```

Convert linear magnitude to decibels: 20*log10(|x|).
