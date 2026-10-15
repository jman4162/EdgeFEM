# Changelog

## [Unreleased]

### Added

- **Optional SuiteSparse UMFPACK backend for direct solves.** CMake option
  `EDGEFEM_USE_SUITESPARSE` (default OFF). With `-DEDGEFEM_USE_SUITESPARSE=ON`
  the build links UMFPACK found by `cmake/FindUMFPACK.cmake`
  (`brew install suite-sparse`; Debian/Ubuntu `libsuitesparse-dev`) and
  configure fails if it is missing. `SolveOptions.direct_backend`
  (`DirectBackend.Default`, `SparseLU`, `UmfPack`) picks one per solve;
  `umfpack_available()` reports the build. UMFPACK runs with its
  METIS_GUARD ordering. `DirectLU` (C++) is the factorization behind
  `solve_linear`, `solve_linear_multi`, `calculate_sparams*` and
  `frequency_sweep`.

  Licensing: UMFPACK is licensed GPL-2.0-or-later. An ON build links
  `libumfpack` (and through it CHOLMOD, AMD, the other SuiteSparse
  libraries and libomp), so a distributed binary or wheel from an ON build
  contains or depends on GPL-licensed code and carries the GPL's
  conditions; EdgeFEM's own source stays MIT. The default (OFF) build links
  no SuiteSparse library and uses Eigen SparseLU. On the 10 GHz patch,
  SparseLU factors in 2.2 s (density 6) and 22 s (density 12), UMFPACK in
  0.19 s and 0.73 s.
- **`SweepResult.min_reflection_index`** (C++ `std::optional`, Python `None`
  or int): with `keep_min_reflection`, the frequency whose fields were kept.
  `StackedPatchDesign.frequency_sweep()` and `PatchAntennaDesign`'s FEM sweep
  return a `FrequencySweepResult`, which unpacks as
  `(frequencies, Z_in, S11)` as before and adds `min_s11_index` and
  `min_s11_frequency`.
- **`FreeDofs`** (`edgefem/free_dofs.hpp`): numbering of the non-PEC edges,
  with `restrict_matrix`, `restrict_vector` and `expand`.
- **`tests/bench_direct_lu`** (built with `EDGEFEM_BUILD_DEBUG_TOOLS`):
  factor + solve time and peak memory for one backend and ordering.

- **`frequency_sweep(..., keep_fields=...)`** on `StackedPatchDesign` and
  `PatchAntennaDesign` keeps field solutions from the sweep's own
  factorizations: `"min_s11"` (the field at the smallest |S11|; the sweep
  holds only the best so far), a list of sweep frequencies in Hz, or both.
  `radiation_pattern()`, `directivity()`, `near_field()` and
  `export_fields_vtk()` at a kept frequency then run no solve. Default:
  nothing kept, as before. The C++ `frequency_sweep` takes a
  `SweepFieldRequest` and returns `SweepResult::fields`; Python:
  `keep_indices`, `keep_min_reflection`.
- **`calculate_sparams_with_fields`** (C++ and Python): the S-matrix and the
  per-port field solutions from one factorization.
- **`Mesh.node_coords()`, `tet_nodes()`, `tet_edges()`, `tet_edge_orient()`,
  `tet_phys()`, `tri_nodes()`, `tri_phys()`**: numpy copies of the mesh,
  built once per call. The list attributes (`nodes`, `tets`, `tris`) still
  build a Python object per element on every access.

### Changed

- **Direct solves drop the PEC rows and columns.** `calculate_sparams*` and
  `frequency_sweep` factor the system restricted to the non-PEC edges and
  expand each field to one entry per mesh edge (zero on PEC edges); the
  full-length fields equal the former ones to 1e-9 relative
  (`tests/test_direct_solvers.cpp`). The identity rows kept their explicit
  zeros in the sparsity pattern, so the orderings saw a denser matrix.
  10 GHz RO4350B patch, one factorization (`tests/bench_direct_lu`):

  | | density 6 (31,924 edges) | density 12 (71,454 edges) |
  |---|---|---|
  | before: SparseLU, full system | 10.8 s, 2.2 GB | 58.1 s, 6.3 GB |
  | SparseLU, PEC eliminated | 2.2 s, 0.8 GB | 22.2 s, 3.7 GB |
  | UMFPACK, PEC eliminated | 0.19 s, 0.24 GB | 0.73 s, 1.0 GB |

  The SparseLU fallback keeps COLAMD: AMD and METIS orderings took 478 s
  and 282 s at density 6.
- **Factorization counter renamed** to `direct_factorization_count()`
  (Python `_direct_factorization_count`); it counts every backend.
- **`StackedPatchDesign.simulate()` factors the system once** (was twice:
  `calculate_sparams` and a second solve for the field). At 10 GHz on the
  RO4350B patch at density 6 (25,260 tets) `simulate()` takes 10.4 s
  instead of 22.4 s; the field equals the former one to 1e-10 relative.
- **`near_field()`** locates points with numpy barycentric coordinates over
  all tets and confirms the first candidate with the C++ test, instead of
  copying the node list four times per tet visited. Results are unchanged.
- **Pattern cache bounds.** Storing a new field at a frequency drops that
  frequency's cached patterns, the cache holds at most 8 patterns (least
  recently used dropped first), and a cached pattern is used only while
  the field it came from is still the stored one. `_last_solutions` is a
  read-only view; `_store_solution()` is the only write path.
- **`StackedPatchDesign.radiation_pattern` is vectorized and cached.** The
  Stratton-Chu sum over Huygens faces runs as blocked numpy matrix products
  instead of a Python loop over directions and faces (results agree with the
  loop to 1e-15 relative to the peak). Patterns are cached per
  (frequency, n_theta, n_phi) and dropped when the field solution at that
  frequency is replaced or the mesh is regenerated. `directivity()` takes
  optional `n_theta`/`n_phi` (default 91 x 73, as before) and reuses a
  cached pattern when the grid matches.

### Fixed

- **A failed solve returns no field.** `calculate_sparams_with_fields` and
  `frequency_sweep` returned the last iterate of an iterative solve that
  did not converge as its field, while `SParamSolution::x` and
  `SweepField::x` were documented as empty in that case. They are empty
  now, the S-matrix column stays NaN, and `StackedPatchDesign.simulate()`
  stores no field for a failed solve (`radiation_pattern()` and the other
  field users raise `RuntimeError` instead of reading a wrong field).
- **Stratton-Chu far-field kernel: three defects, one fix.** (1) The
  radiation-integral kernel applied `jk0` to the magnetic-current term only;
  both current terms carry it, so J was mis-weighted against M by roughly a
  factor of k0 and patterns from a Huygens surface were effectively M-only.
  (2) The J term used `-Z0 rhat x J` where the kernel requires
  `Z0 (rhat x J) x rhat`, and the M term used the transverse projection where
  it requires `rhat x M`. (3) Eigen's `cross()` conjugates complex operands'
  products, silently corrupting every complex cross in the chain — the file
  now uses an explicit conjugation-free cross, and pattern projection no
  longer conjugates the field (Eigen `dot()` conjugates its first argument).
  The kernel also gains the physical `1/(4 pi)` prefactor, making the output
  an absolute r-normalized far field.

  **Pattern outputs shift.** Ratio-based metrics (directivity, HPBW,
  normalized pattern shape) were approximately right before because the
  dominant M term alone carries most of the shape for outgoing fields;
  cross-polarization and any absolute quantity (RCS, gain from absolute
  fields) were wrong. `tests/test_ntf.cpp` now checks the kernel against the
  analytic Hertzian-dipole far field — absolute magnitude to 2%, pure
  polarization, and directivity 1.5 — with assertions that survive
  `-DNDEBUG` builds (the previous asserts silently vanished in release
  builds). Prerequisite for the planned RCS capability.

## Unreleased (2026-07-14) — assembled port operator, mesh-convergent S-parameters

Wave-port formulation rework in `calculate_sparams_eigenmode`; no API
breaks (one new builder, one changed default).

### Changed

- Port ABC is now the assembled surface mass matrix: `A += jβ·M_s` with
  M_s = ∫(n̂×N_i)·(n̂×N_j)dS over the port triangles, replacing the
  diagonal-lumped operator. Excitation, normalization, and modal
  extraction all use the M_s inner product.
- Port modes come from a 2D generalized eigensolve (K_s, M_s) on the port
  face (`build_wave_port_2d` / `solve_port_mode_2d`), replacing the 3D
  cavity eigensolve. The discrete cutoff replaces the analytical one so β
  is consistent with the profile. Port setup at 12k tets: 1014 s → <0.1 s.
- `port_abc_scale` default 0.5 → 1.0. The sweep optimum is measured at
  exactly 1.0; the knob is diagnostic only. Per-port β replaces the
  port-0-only β (latent multi-port bug), and β is now computed from the
  material adjacent to the port (complex for lossy fills), so
  dielectric-filled guides get a matched port impedance.

### Fixed

- Triangle local-edge-table mismatch: the mesh loader fills Tri3 edges in
  order (0,1),(1,2),(2,0), but `lumped_port.cpp` (and initially the new
  surface assembly) used (0,1),(0,2),(1,2), silently scattering surface
  integrals to wrong edges. All surface code now uses the loader's order;
  `test_triangle_mass_matrix` guards the convention.

### Measured impact (WR-90)

- Band sweep 7-12 GHz: |S11| 0.003-0.053 (was 0.04-0.14), |S21| ≥ 0.997.
- Mesh convergence at 10 GHz: |S11| 0.15 → 5×10⁻⁴ and |S21| error
  1.3% → 0.01% from 5 to 18 elements/wavelength. The previous 1-2%
  magnitude floor (non-convergent) is eliminated; phase error unchanged
  at ~O(h²).
- Full CTest suite: 347 s → 62 s.

## Unreleased (2026-07-11) — documentation and validation integrity pass

Outcome of a full technical review of v1.0.0. No solver physics changed;
this pass makes the documentation, tests, and examples describe what the
code actually does.

### Documentation

- README, docs, and CLAUDE.md now state implemented capability: single-axis
  Bloch periodic BCs (not "Floquet ports"), analytically validated band of
  roughly 5–26 GHz (not "100 kHz–110 GHz"), single-threaded execution (no
  HPC claims), measured accuracy ranges (no single-point "99.4%" framing).
- Current Limitations expanded: the "PML" is a graded absorber using a
  scalar average of the stretch tensor; wave ports are first-order lumped
  ABCs with a tuned 0.5 scale factor; no TEM/inhomogeneous ports; periodic
  constraint elimination is dense; Gmsh v2 ASCII input only.
- `system_design_spec.md` and `requirements.md` marked as aspirational.
- Removed `SPEC_TIME_DOMAIN.md` (unedited chat transcript),
  `project_review.md` (stale generated review), and
  `docs/VALIDATION_REPORT.md` (conflicted with `docs/validation.md` on the
  same benchmark). Fixed placeholder citation in `papers/README.md`;
  rewrote `AGENTS.md` to list only tools that exist.

### Validation

- `docs/validation.md` rebuilt from measured data only, with a reproduce
  command per table (WR-90 7–12 GHz, ABC scaling sweep, mesh convergence,
  cavity eigenmodes, lossy attenuation) and an explicit "What is NOT
  validated" section.
- `scripts/run_convergence_study.py` now runs the full solve chain and fits
  the observed convergence order (the previous version generated meshes but
  returned `None` from `run_simulation()`). Measured WR-90 result: phase
  error converges ~O(h²) (18.6° → 0.9°); S-magnitude error plateaus at a
  1–2% floor set by the lumped port ABC and does not converge with mesh.
- Fresnel slab test computes R+T = 1 instead of hard-coding the result;
  it and the coax test are labeled setup-sanity tests (neither compares an
  FEM solve against the analytic reference).
- 22 assertion-free diagnostic binaries moved behind
  `EDGEFEM_BUILD_DEBUG_TOOLS` (default OFF).
- Fixed `VECTOREM_PYTHON` → `EDGEFEM_PYTHON` in tests/CMakeLists.txt: the
  `python_smoke` CTest had never been registered, in CI or locally. It now
  runs with the interpreter the bindings were built against.

### Fixed

- Bound `compute_barycentric` in pyedgefem; `StackedPatchDesign.near_field()`
  previously raised `AttributeError`.
- `tools/agents/mesh_qa.py` is a real mesh checker (tet volumes, aspect
  ratios, Gmsh format detection) instead of a no-op that always reported ok.
  It caught and removed an unreadable Gmsh v4 mesh shipped in examples/.
- `examples/run_patch_fullwave.py` runs the actual FEM patch sweep
  (previously assembled with no BCs and no ports).
- `examples/unit_cell_demo.py` runs the real periodic solve (previously
  printed a fabricated reflection trend).
- `gmsh` declared as a runtime dependency in pyproject and CI requirements.
- `UnitCellDesign.load_mesh` double-load; dead methods and a stale warning
  in `PatchAntennaDesign`.

### Removed from version control

- `.DS_Store`, generated Touchstone/CSV/mesh outputs, and the `simulations/`
  scratch directory (~75k lines of committed artifacts).

## 1.0.0 (2026-02)

Initial release: lowest-order Nédélec edge-element frequency-domain solver
with wave/lumped ports, single-axis periodic BCs, dispersive material
models, Stratton-Chu far fields, Python SDK, and PyPI packaging.
