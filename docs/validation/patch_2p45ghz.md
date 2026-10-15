# Probe-fed patch: mesh convergence and openEMS cross-check

Measured 2026-10-15 on macOS ARM64 (14 cores, 36 GB), Release build with
UMFPACK (`EDGEFEM_USE_SUITESPARSE=ON`). EdgeFEM at commit 13cab87 (surface
mass ABC) with the Python defaults of 23c6271; openEMS from the
`antenna-cad-openems` Docker image (openEMS-Project 84703e2).

## Geometry

The 2.45 GHz FR-4 patch of the known-answer test
(`python/tests/test_patch_physics.py`), built by `PatchAntennaDesign`:

| Item | Value |
|---|---|
| Patch | L 28.6 mm (y), W 37.3 mm (x), PEC sheet at z = h |
| Substrate | FR-4, h 1.6 mm, eps_r 4.4, tan d 0.02, 111.9 x 111.9 mm |
| Ground | PEC at z = 0, 111.9 x 111.9 mm (3 x W), same size as the substrate |
| Probe | radius 0.5 mm, at y = -9 mm from the patch center: square PEC post of half-side 0.424 mm (equivalent radius 0.59 x side) |
| Ground hole | square, half-side 0.975 mm (hole/probe radius ratio 2.30) |
| Port | 50 ohm gap port on the strip from the post face to the hole edge at z = 0, E along +x |
| Cavity model (Balanis 14.2, Hammerstad dL) | 2.4669 GHz |

## EdgeFEM convergence

Reproduce: `PYTHONPATH=build-fix/python:python python tools/patch_convergence.py --case fr4 --density 8 12 16`
and `... --density 12 --patch-factor 4 3`.
Air box lambda/2 at 1.8 GHz. f_res is the |S11| minimum of a 2.25 to
2.60 GHz sweep in 10 MHz steps, refined by a parabola through the three
points around it; Re(Z) off resonance comes from a 1.8 to 3.0 GHz sweep
in 100 MHz steps. The near-patch size is
`_PATCH_SIZE_X_DENSITY * h / density`. Time is the two sweeps (49 LU
factorizations) plus meshing; memory is the peak resident set of the
process.

| Density | Near-patch size | Edges | f_res (GHz) | vs cavity | min \|S11\| (dB) | Z at f_res (ohm) | Re Z 1.8 / 3.0 GHz (ohm) | Time | Peak |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.75 h | 63,950 | 2.3125 | -6.26% | -15.24 | 63.8 - j14.5 | 0.99 / 1.66 | 32 s | 7.5 GB |
| 12 | h/2 | 121,532 | 2.3600 | -4.34% | -14.53 | 60.3 - j18.3 | 0.50 / 1.40 | 97 s | 7.2 GB |
| 16 | 0.375 h | 225,156 | 2.3703 | -3.92% | -14.39 | 60.8 - j18.5 | 0.47 / 1.35 | 295 s | 8.9 GB |
| 12 | h/3 (default) | 258,355 | 2.3768 | -3.65% | -14.39 | 55.7 - j19.9 | 0.45 / 1.34 | 303 s | 9.7 GB |
| 12 | h/4 | 456,162 | none | | | | | | |

At h/4 every solve failed: UMFPACK reported "factorization failed" with
4.2 GB resident (SUSPECTED: the 32-bit-index UMFPACK interface Eigen uses
for `SparseMatrix<complex<double>>` runs out of addressable workspace).

f_res rises as the mesh near the patch shrinks: +47.5 MHz from 0.75 h to
h/2, +16.8 MHz from h/2 to h/3 (step ratio 1.5 both times). The ratio of
the two increments, 2.83, gives an observed order of 2.6 and an
extrapolated f_res of 2.386 GHz (SUSPECTED: three points, one of them on
the coarse side). The default near-patch size is h/3 at density 12
(`_PATCH_SIZE_X_DENSITY = 4.0`), 0.4% below that extrapolation; h/2 was
1.1% below at a third of the cost.

The known-answer test (25 MHz steps, air box lambda/2 at 2.25 GHz) finds
2.375 GHz, -3.73% from the cavity model, min |S11| -14.39 dB, in 112 s and
8.4 GB (`pytest python/tests -m slow`). It found 2.350 GHz (-4.74%) with
the h/2 default.

## ABC: surface mass matrix

Commit 13cab87 replaced `j k0` on the diagonal of each ABC edge with
`j k0` times the assembled surface mass matrix. Density 12, near-patch
h/2, same sweeps as above:

| Air box | ABC | f_res (GHz) | min \|S11\| (dB) | Re Z 1.8 GHz | Re Z 3.0 GHz |
|---|---|---|---|---|---|
| lambda/2 at 1.8 GHz (83 mm) | diagonal | 2.3586 | -15.11 | 0.49 | 1.38 |
| lambda/2 at 1.8 GHz (83 mm) | surface mass | 2.3600 | -14.53 | 0.50 | 1.40 |
| 30.6 mm (lambda/4 at 2.45 GHz) | diagonal | 2.3652 | -12.70 | 0.49 | 1.27 |
| 30.6 mm (lambda/4 at 2.45 GHz) | surface mass | 2.3621 | -14.55 | 0.52 | 1.36 |

CONFIRMED (table): with the surface mass ABC, the lambda/4 and lambda/2
boxes agree to 2.1 MHz in f_res and 0.02 dB in min |S11|; with the
diagonal ABC they differed by 6.6 MHz and 2.4 dB. The off-resonance
Re(Z) of 16.5 ohm (1.8 GHz) and 26.4 ohm (3.0 GHz) in the 2026-10-12
investigation came from the model before the gap-port and geometry fixes;
on the current model it is 0.5 and 1.4 ohm with either ABC, and openEMS
gives 0.38 to 0.53 and 1.08 to 1.20 ohm (gap-port runs below).

## openEMS cross-check

Reproduce: `python3 tools/openems_patch_crosscheck.py --case fr4 --workdir <dir>`
(runs openEMS in Docker; see the script header). Same patch, substrate,
ground, post, hole and gap port as EdgeFEM; 1.8 to 3.0 GHz Gaussian
excitation, 1 MHz output steps, PML_8 on all six sides, lambda/4 of air
around the board, end criterion -50 dB (reached in every run), thirds rule
at the patch edges.

| Run | Patch-edge step | Cells across h | Cells | f_res (GHz) | min \|S11\| (dB) | Z at f_res (ohm) | f at max Re Z (GHz) | Time |
|---|---|---|---|---|---|---|---|---|
| finite ground | 0.4 mm | 4 | 435,891 | 2.379 | -12.95 | 68.1 - j20.0 | 2.356 | 285 s |
| finite ground | 0.2 mm | 8 | 576,300 | 2.389 | -13.03 | 66.6 - j20.6 | 2.364 | 623 s |
| ground and substrate through the PML | 0.4 mm | 4 | 427,893 | 2.379 | -13.17 | 67.8 - j19.2 | 2.356 | 246 s |
| vertical lumped port, no post or hole | 0.4 mm | 4 | 435,891 | 2.392 | -16.89 | 56.6 - j13.9 | 2.362 | 261 s |

EdgeFEM at the default mesh (density 12, h/3) against openEMS on the finer
mesh:

| | EdgeFEM | openEMS | Difference |
|---|---|---|---|
| f_res | 2.3768 GHz | 2.389 GHz | -0.51% |
| min \|S11\| | -14.39 dB | -13.03 dB | 1.4 dB deeper |
| Z at f_res | 55.7 - j19.9 ohm | 66.6 - j20.6 ohm | |
| f at max Re Z | 2.35 GHz (10 MHz grid) | 2.364 GHz | -0.6% |

CONFIRMED (tables): f_res agrees to 0.51%, inside the 2% target. At
density 12 with h/2 near the patch EdgeFEM is 1.2% low (2.3600 GHz). Both
solvers put the resonance 3.2% (openEMS) to 3.65% (EdgeFEM) below the
cavity model; the vertical-port feed moves openEMS by +0.55% in f_res and
+0.25% in the Re Z peak, so the gap from the cavity model is not a feed
artifact (SUSPECTED cause: the
Hammerstad length extension and eps_reff for this W/h).

Where the two models differ:

- Domain truncation. EdgeFEM's domain is the ground-plane footprint: the
  ABC side walls stand on the ground edges and there is nothing below the
  ground. openEMS has a finite ground in air with lambda/4 of air on all
  sides. Extending the openEMS ground and substrate through the PML
  changed f_res by 0 MHz and min |S11| by 0.2 dB, so the ground edge does
  not drive the difference here.
- The ground hole outside the port strip is a natural (PMC) boundary in
  EdgeFEM and open to the air below in openEMS.
- Loss: EdgeFEM uses eps_r (1 - j tan d) at every frequency; openEMS a
  conductivity kappa = 2 pi f0 eps0 eps_r tan d with f0 = 2.45 GHz, so
  tan d is 0.02 only at 2.45 GHz (0.0206 at 2.38 GHz).
- Discretization: tetrahedral first-order edge elements versus a
  staircase-free Cartesian FDTD grid (every face is axis-aligned).

## EdgeArray's 10 GHz RO4350B patch

`tools/patch_convergence.py --case ro4350b_10g`: L 7.679 mm, W 9.820 mm,
h 0.508 mm, eps_r 3.66, tan d 0.0037, ground 24.81 mm, probe radius 0.5 mm
at the inset point (y = -L/2 + 2.131 mm), as EdgeArray's EdgeFEM worker
builds it; `f_min` 8.5 GHz. Cavity model 10.000 GHz. The sweep is
EdgeArray's live-test sweep, 8.5 to 11.5 GHz in 100 MHz steps.

| Density | Near-patch | Edges | min \|S11\| (dB) at | Interior argmin | Phase spread | f_res vs 10 GHz | Time |
|---|---|---|---|---|---|---|---|
| 6 | h/1.5 | 45,124 | -13.77 at 9.0 GHz | yes | 364 deg | -10.0% | 34 s |
| 10 | 0.4 h | 97,114 | -14.09 at 9.3 GHz | yes | 367 deg | -7.0% | 121 s |
| 12 | h/3 | 145,912 | -19.31 at 9.3 GHz | yes | 368 deg | -7.0% | 278 s |

With 25 MHz steps, f_res is 9.046, 9.254 and 9.296 GHz at density 6, 10
and 12. openEMS on the same model gives 9.510 GHz (-17.35 dB); with a
0.25 mm probe 9.593 GHz, with a vertical lumped port 9.620 GHz. So the
probe-fed model of this patch resonates 4 to 5% below the cavity model in
openEMS, and EdgeFEM at density 12 is another 2.3% below openEMS. The
case meets min |S11| < -6 dB, interior argmin and phase spread > 30 deg;
it fails f_res within 5% of 10 GHz at every density, and would still fail
if EdgeFEM matched openEMS exactly (9.51 GHz is -4.9%, on the edge).
