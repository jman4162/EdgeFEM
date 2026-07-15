# Changelog

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
