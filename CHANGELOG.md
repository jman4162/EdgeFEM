# Changelog

## 1.1.0 (2026-10-15)

### Breaking changes and known limitations

- **Modal and Floquet port S-parameters are not validated (known
  limitation).** Only ports built by `build_lumped_port` get the new
  gap-port terms; `WavePort.kind` (`PortKind.Lumped` / `PortKind.Modal`,
  default Modal) selects them. Ports from `build_wave_port`,
  `build_wave_port_2d`, `build_wave_port_from_eigenvector` and hand-built
  `WavePort`s keep the 1.0.0 terms (load w w^H / Z0, source
  2 w / sqrt(Z0), S(j, i) = V_j / sqrt(Z0_i) - delta_ij) in
  `assemble_maxwell`, `calculate_sparams(_with_fields)`, `frequency_sweep`
  and `calculate_sparams_periodic`. On the same meshes they equal 1.0.0 to
  1e-13 (`test_modal_port_1p0`, and `test_unit_cell_design.py` to 1e-9).
  Those 1.0.0 values are not physical: WR-90 with analytic TE10 weights
  gives |S11| = 2.41 / 1.63 / 1.33 at 8 / 10 / 12 GHz, and a bare grounded
  slab (eps_r 3.5, 1 mm, port plane 14.99 mm above) gives arg(R) =
  -172.15 / -175.12 / -176.70 deg at 8 / 10 / 12 GHz against the analytic
  -127.67 / +155.07 / +77.58 deg. These functions print a one-time stderr
  warning (C++) and emit a `RuntimeWarning` (Python) when a Modal port is
  present. `calculate_sparams_eigenmode` is unchanged, validated on WR-90
  and WR-42, and does not warn. A hand-built `WavePort` meant as a gap port
  needs `kind = PortKind.Lumped`.
- **`assemble_maxwell` with a lumped port.** The active port's right-hand
  side is now the gap-port source b = j k0 eta0 (2 / sqrt(Z0)) v, and every
  lumped port adds j k0 eta0 / Z0 v v^H to A. In 1.0.0 it was
  b = 2 w / sqrt(Z0) and A += w w^H / Z0 (with build_lumped_port scaling
  ||w|| to sqrt(Z0)). For one edge with unit weight and Z0 = 50 ohm, b goes
  from 0.2828 to j 2.233e4 at 10 GHz (j 3.5543e-7 at omega = 1 rad/s), so a
  caller that solves the assembled system itself gets a different field
  scale and phase than in 1.0.0. The solution now corresponds to an
  incident power wave a = 1 and the port voltage is V = v^H x.
- **`build_lumped_port` weights** are the voltage functional v = w / W
  (W = strip width) instead of w scaled to ||w|| = sqrt(Z0), and it throws
  (`RuntimeError` in Python) when `e_direction` has no extent in the port
  surface.
- **`normalize_port_weights` is a no-op** for every port kind. It prints a
  one-time stderr warning (C++) or raises a `DeprecationWarning` (Python)
  and returns the ports unchanged. In 1.0.0 it rescaled the weights to
  w^H A^-1 w = Z0 (in place in C++, as a returned copy in Python). It will
  be removed in a later release.
- **`StackedPatchDesign` raises on geometry it cannot build.** `ValueError`
  when `probe_hole_radius <= probe_radius` (constructor), and from
  `generate_mesh()` when the probe does not fit inside the first patch or
  the probe hole does not fit inside the ground plane. `generate_mesh()`
  raises `RuntimeError` for a surface it cannot classify (ground, port,
  probe, patch, ABC, Huygens), for a port strip whose area differs from
  (b - a) x 2a by more than 1e-6 relative, for a group with no surfaces,
  and for a substrate layer with no volume.
- **`generate_mesh(f_min=None)` default.** The air box is lambda/2 at
  `f_min`, and `f_min` defaults to 0.8 x `design_freq`. 1.0.0 used lambda/4
  at the design frequency, so the default box is 2.5 times taller and the
  mesh larger.
- **`cavity_size` and `cavity_depth` are deprecated and ignored** (a
  `DeprecationWarning`); the model has no cavity below the ground plane and
  the cavity physical tags are gone.
- **Other changes visible to direct callers.**
  - `use_abc` adds j k0 times the boundary surface mass matrix instead of
    j k0 on the diagonal of each ABC edge, for every port kind
    (see Changed).
  - A failed solve returns an empty field and a NaN S-matrix column;
    `StackedPatchDesign.frequency_sweep()` gives NaN for `Z_in` there (was
    1e6) and a `RuntimeWarning`.
  - Radiation patterns from a Huygens surface change with the Stratton-Chu
    fix (see Fixed): cross-polarization and absolute levels differ from
    1.0.0.
  - C++: `frequency_sweep` takes a `SweepFieldRequest` (defaulted) and
    `check_solver_convergence` moved to `solver.hpp`; code built against
    1.0.0 headers must be recompiled.
  - The test-only factorization counter is `direct_factorization_count()`
    (Python `_direct_factorization_count`).

### Added

- **openEMS cross-check and patch convergence study.**
  `tools/openems_patch_crosscheck.py` builds the probe-fed patch of
  `StackedPatchDesign` (post, ground hole, gap port, finite ground) in
  openEMS and runs it in the `antenna-cad-openems` Docker image;
  `docs/validation/patch_2p45ghz.md` has the method and numbers. 2.45 GHz
  FR-4 patch: EdgeFEM 2.3768 GHz (density 12), openEMS 2.389 GHz, cavity
  model 2.4669 GHz.
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

- **UMFPACK runs with 64-bit indices** (`umfpack_zl_*`; `DirectLU` hands it
  a `SparseMatrix<complex<double>, ColMajor, SuiteSparse_long>` copy). With
  32-bit indices the 2.45 GHz FR-4 patch at density 12 with h/4 near the
  patch (456k edges) failed with `UMFPACK_ERROR_out_of_memory` on a 76 GB
  symbolic estimate; with 64-bit indices it factors in 20 s and 6.4 GB and
  gives f_res 2.388 GHz (openEMS 2.389 GHz). Peak memory of other runs rises
  by 1.4 to 2.2 GB (known-answer test 8.4 -> 10.6 GB, same S11).
- **Patch mesh: finer near the patches, larger default air box.**
  `StackedPatchDesign._PATCH_SIZE_X_DENSITY` 6.0 -> 4.0: the mesh size
  inside and near the patch footprint is h/3 at density 12 (was h/2). On
  the 2.45 GHz FR-4 patch at density 12 f_res moves 2.360 -> 2.377 GHz
  (openEMS: 2.389 GHz) and the solve time triples (97 s -> 303 s, peak
  7.2 -> 9.7 GB); the known-answer test now finds 2.375 GHz (-3.7% from the
  cavity model, was -4.7%) in 112 s. `generate_mesh(f_min=None)` now sizes
  the air box for 0.8 x design_freq (was design_freq), so a sweep within
  +/-20% of design_freq raises no air box warning; the fast Python suite
  went from 30 warnings to 13 (all from the analytical
  `PatchAntennaDesign` paths) and from 68 s to 123 s.
  `tools/patch_convergence.py` reproduces the study.
- **First-order ABC uses the boundary surface mass matrix.** `use_abc` now
  adds j k0 int_S N_i . N_j dS over the ABC triangles (the weak form of
  n x curl E + j k0 n x n x E = 0) instead of j k0 on the diagonal of each
  ABC edge, in `assemble_maxwell` and `frequency_sweep` alike.
  `abc_surface_matrix()` returns the frequency-independent matrix and
  `assemble_surface_mass()` generalizes `assemble_port_surface_mass()` to a
  set of tags. `test_abc` checks the assembled term entry by entry and is no
  longer `known_issue`. On the 2.45 GHz FR-4 patch (density 12, 10 MHz
  steps) the air-box height mattered less afterwards: with the box lambda/2
  at 1.8 GHz, f_res 2.3586 -> 2.3600 GHz and min |S11| -15.11 -> -14.53 dB;
  with a lambda/4 box, 2.3652 -> 2.3621 GHz and -12.70 -> -14.55 dB.
  Off-resonance Re(Z) was already small after the gap-port and geometry
  fixes and did not change: 0.49 -> 0.50 ohm at 1.8 GHz, 1.38 -> 1.40 ohm
  at 3.0 GHz.
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
  cached pattern when the grid matches. 10 GHz RO4350B patch, 46 x 73
  directions: the Python loop took 14.2 s (90 Huygens faces, density 6) and
  33.7 s (162 faces, density 12); the vectorized sum takes 0.01 s and 0.02 s,
  with E_theta and E_phi equal to the loop to about 1e-15 relative to the
  peak.

### Fixed

- **The lumped port is a gap port with a physical load and source.**
  The old port scaled its weights to ||w|| = sqrt(Z0) and added a
  frequency-independent real load, without the j k0 eta0 factor that
  relates a gap current to the field equation, so every S11 sat near -1.
  The port now stores v = w / W (W = strip width),
  so V = v^H x is the gap voltage; the load is the term
  j k0 eta0 / Z0 v v^H added to the system matrix; the source is a Norton
  current 2 / sqrt(Z0) in parallel with Z0 (incident power wave a = 1).
  S(j, i) = V_j / sqrt(Z0_j) - delta_ij, and all ports share one loaded
  matrix, so `calculate_sparams` and `frequency_sweep` factor once per
  frequency (was once per port) and solve the right-hand sides together
  (`solve_linear_multi`). Derivation: comment block at the top of
  `src/assemble_maxwell.cpp`. `normalize_port_weights` is a no-op that warns
  (`DeprecationWarning` in Python); `StackedPatchDesign` no longer calls it.

  Before and after (2 / 5 / 8 GHz unless noted):
  - Closed lossless box: S11 was -0.9938 / -0.9936 / -0.9933, a short that
    ignores the box; it is now -0.986+0.168j, -0.910+0.415j, -0.762+0.647j,
    |S11| = 1, Z = j4.23, j10.87, j18.36 ohm.
  - TEM line: |S11| was 1.07 / 0.98 / 1.09 and |S21| 0.11 / 0.03 / 0.09; it
    is now |S11| = 6.1e-5, 4.9e-4, 6.5e-5 and |S21| = 1, with Z11 within
    2.3e-5, 1.6e-3, 9.1e-3 of -j Zc cot(k0 L).
  - WR-90 with analytic TE10 weights through `calculate_sparams`:
    |S11|^2 + |S21|^2 was 12.3 / 5.8 / 2.65 at 7 / 8 / 10 GHz; with the gap
    terms applied to these modal ports it was 1.0000 from 7 to 100 GHz.
    Modal ports no longer get the gap terms (see Breaking changes), so these
    sums are back to their 1.0.0 values.
  - Probe-fed patch (density 6, 2.2 / 2.45 / 2.7 GHz): S11 was -0.9947 /
    -0.9938 / -0.9921; it is -0.708+0.686j, -0.560+0.814j, -0.288+0.934j.
  - `test_sparam_accuracy` takes 45.5 s (was 72.5 s) and `test_frequency_sweep`
    2.2 s (was 6.6 s).

  New tests: `test_lumped_lossless_box`, `test_lumped_tem_line`; the
  `test_wave_port_excitation` expected values changed to the physical source
  (j 3.5543e-7 at omega = 1 rad/s, Z0 = 50 ohm).
- **Patch geometry: substrate volumes, no cavity, probe from
  `probe_radius`.** `StackedPatchDesign` built no substrate volume
  (the substrate tag held only the probe interior, 41.6 mm^3 against
  20,034.6 mm^3 expected for the FR-4 patch, so the substrate was solved as
  air), and the PEC probe ran down to a PEC cavity floor at z = -5 mm,
  shorting it under the port. Each substrate layer is now its own
  volume (physical tag 110 + layer index), the domain starts at the ground
  plane (z = 0), and the lumped port is a rectangle in a square hole in the
  ground, from the probe's +x face to the hole edge. The probe is a square
  PEC post of side `probe_radius / 0.59` (equivalent thin-wire radius
  `probe_radius`) at every mesh density. New keyword `probe_hole_radius`
  (default 2.30 x `probe_radius`, the 50 ohm air-coax ratio). `cavity_size`
  and `cavity_depth` are ignored and raise a `DeprecationWarning`; the
  cavity tags are gone. The air box is `0.5 c0 / f_min` above the top metal
  (new `generate_mesh(f_min=...)`, new property `air_box_height`), and
  `simulate()` / `frequency_sweep()` warn once per mesh when run below the
  frequency at which the box is lambda/2.

  Checked on the 2.45 GHz FR-4 patch at density 12: the substrate volume is
  20,033.43 mm^3 against 20,034.58 mm^3 nominal (the difference is the
  1.15 mm^3 probe post), and the ground, patch, probe-wall and port areas
  equal their nominal values to 1e-6. The known-answer test (cavity model
  2.4669 GHz, 5% tolerance) previously saw |S11| of -0.02 to -0.04 dB over
  2.25 to 2.70 GHz with its minimum at the sweep edge; it now finds an
  interior minimum at 2.350 GHz (-4.74%, -13.66 dB) with the geometry as
  first committed, and 2.375 GHz (-3.73%) with the mesh default below.
  A closed lossless box with the patch gives |S11| = 1.000000 at 13 points
  from 2.30 to 2.45 GHz.
- **A failed sweep point no longer becomes a number.**
  `StackedPatchDesign.frequency_sweep()` mapped a NaN S11 (failed solve) to
  `Z_in = 1e6`; it now returns NaN for Z_in there and raises a
  `RuntimeWarning` naming the failed frequencies. The C++
  `frequency_sweep` warns on stderr for each failed solve, as
  `calculate_sparams` does (`check_solver_convergence`, now in
  `solver.hpp`), and `DirectLU` failure messages carry the backend's status:
  `UMFPACK factorization failed, status -1 (UMFPACK_ERROR_out_of_memory);
  symbolic estimates: ...` or SparseLU's own message.
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

### Build and packaging

- **Version 1.1.0, one source.** `__version__` in `python/edgefem/__init__.py`
  is the only place the number is written; `pyproject.toml` takes it through
  scikit-build-core's regex metadata provider and `CMakeLists.txt` reads the
  same line for `project(... VERSION ...)`.
- **`uv build --wheel` / `pip wheel .` work again.** `pyproject.toml` used
  `cmake.minimum-version` and a table-form `wheel.packages`, which
  scikit-build-core 0.8 and later reject ("Use cmake.version instead of
  cmake.minimum-version"). It now uses `cmake.version = ">=3.20"` and a list
  for `wheel.packages`, and requires `scikit-build-core>=0.10` (built with
  0.10.7 and 1.1.1) and `pybind11>=2.12` (2.11 predates Python 3.13). Added
  the Python 3.13 classifier.
- **`.github/workflows/wheels.yml`:** cibuildwheel 3.4.1 (was 2.21.3,
  which does not know `CIBW_TEST_SOURCES`; 3.x copies those files into the
  test directory, and its default Linux image is manylinux_2_28),
  `cp313-*` added, `macos-13` replaced by `macos-15-intel` (x86_64) next to
  `macos-14` (arm64), `MACOSX_DEPLOYMENT_TARGET=14.0`, Eigen 3.4.0
  installed in the Linux images (they offer 3.3, and EdgeFEM needs 3.4).
  The wheel test runs `tools/wheel_smoke.py` (version, build flags and a
  WR-90 solve, about 1 s) and `python/tests/test_contract.py` instead of
  only importing `pyedgefem`. Published wheels are built with `EDGEFEM_USE_SUITESPARSE=OFF`
  because UMFPACK is GPL-2.0-or-later; `pyproject.toml` now sets that option
  to OFF explicitly.
- **`.github/workflows/ci.yml`** runs `ctest -LE known_issue`, so a failing
  test that is not labelled `known_issue` fails the job. `tools/lint.sh`
  (clang-format 19 and clang-tidy, run before ctest) passes: the tree is
  formatted with clang-format 19.1.7.
- **Eigen 3.4.x and 5.x.** `cmake/FindEigen3.cmake` asks for Eigen 3.4,
  then 5, in config mode (each Eigen version file accepts only its own
  major version, so a `3.4...5` range fails on 3.4.0). Homebrew's Eigen
  5.0.1 builds and passes ctest and the fast Python suite.
- **pytest marker `known_issue`** for known-answer tests that fail on an
  open defect; each is also `xfail(strict=True)`.

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
