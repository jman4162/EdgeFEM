# Validation Report

This report contains only numbers measured from code in this repository, with
the command that reproduces each table. Platform for the reported runs:
macOS ARM64, Clang, Release build, BiCGSTAB+ILUT solver.

Scope: validation currently covers hollow-waveguide S-parameters, the port
ABC scaling factor, cavity eigenmodes, and lossy-dielectric attenuation, over
roughly 5-26 GHz. See "What is NOT validated" below.

## WR-90 Rectangular Waveguide (7-12 GHz)

Analytical reference: TE10 with |S11| = 0, |S21| = 1, phase(S21) = -βL
(Pozar). Mesh: 4,227 tets (~10 elements/wavelength), length 50 mm.

Reproduce: `./build/tests/benchmark_wr90`

| Freq (GHz) | \|S11\| | \|S21\| | Phase(S21) | Expected | Phase error |
|-----------|-------|-------|------------|----------|-------------|
| 7.0 | 0.064 | 0.997 | -150.3° | -147.1° | 3.2° |
| 8.0 | 0.093 | 0.994 | 80.9° | 84.8° | 3.9° |
| 9.0 | 0.037 | 0.998 | -15.7° | -10.1° | 5.6° |
| 10.0 | 0.087 | 0.994 | -100.6° | -93.3° | 7.3° |
| 11.0 | 0.042 | 0.997 | 179.8° | -170.3° | 10.1°* |
| 12.0 | 0.143 | 0.987 | 103.8° | 116.6° | 12.8° |

*Phase wrapping near ±180°.

Test pass criteria (deliberately loose for CI stability): |S11| < 0.15,
|S21| > 0.90, phase error < 15°, passivity ≤ 1.05. Do not quote the
single best point as the solver's accuracy; the honest summary is
**|S21| within ~1%, |S11| below ~0.15, and phase error growing to ~13° at
the band edge** at this mesh density.

A WR-42 benchmark (20-26 GHz) runs the same checks: `./build/tests/test_wr42_waveguide`.

## Port ABC Scaling Factor

The wave-port boundary uses a first-order ABC with coefficient α·jβ applied
to port-edge diagonal entries (a lumped approximation of the port surface
mass matrix). α was chosen by sweep, not derived.

Reproduce: `./build/tests/test_abc_scaling_sweep` (writes `abc_scaling_sweep.csv`)

| α | \|S11\| | \|S21\| | Phase err |
|-----|-------|-------|-----------|
| 0.40 | 0.297 | 0.953 | 6.8° |
| 0.45 | 0.188 | 0.980 | 7.1° |
| 0.50 | 0.087 | 0.994 | 7.3° |
| 0.55 | 0.012 | 0.998 | 7.3° |
| 0.60 | 0.092 | 0.994 | 7.3° |
| 1.00 | 0.533 | 0.844 | 5.7° |

Measured optimum on this mesh is α ≈ 0.55; the shipped default is
`port_abc_scale = 0.5`. The sharp sensitivity of |S11| to α is a property of
the lumped-ABC port formulation itself — treat |S11| below ~0.05 from this
port model as formulation-limited, not physical.

## Mesh Convergence (WR-90 at 10 GHz)

Reproduce: `PYTHONPATH=build/python:python python3 scripts/run_convergence_study.py`

Measured (macOS ARM64, July 2026):

| Elements/λ | Tets | \|S21\| | \|S21\| error | \|S11\| | Phase error | Runtime |
|-----------|------|-------|-------------|-------|-----------|---------|
| 5 | 452 | 0.9757 | 2.43% | 0.194 | 18.6° | <0.1 s |
| 8 | 1,323 | 0.9767 | 2.33% | 0.120 | 7.6° | 0.7 s |
| 10 | 2,462 | 0.9943 | 0.57% | 0.094 | 4.0° | 5.3 s |
| 14 | 5,740 | 0.9818 | 1.82% | 0.177 | 2.0° | 95 s |
| 18 | 12,233 | 0.9869 | 1.31% | 0.146 | 0.9° | 1014 s |

Two distinct behaviors, and it is important not to conflate them:

- **Phase error converges cleanly** (18.6° → 0.9°), consistent with the
  expected ~O(h²) dispersion error of lowest-order edge elements. Fitted
  order from these points: ≈2.
- **S-parameter magnitudes do NOT converge with the mesh.** |S21| error and
  |S11| oscillate in the 1-2% / 0.09-0.19 band regardless of density
  (fitted |S21|-error order ≈ 0.45, i.e. no meaningful convergence). The
  magnitude floor is set by the lumped first-order port ABC, not by
  discretization; refining the mesh does not remove it. This is the port
  limitation described in the README, and it is why quoting a single
  best-case |S11| is misleading.

Runtime grows superlinearly (dominated by the dense port eigenvector and
direct-solver fallback); ~10 elements/wavelength is the practical
accuracy/cost point for this formulation.

## Cavity Eigenmodes

Rectangular cavity resonances vs. analytical f_mnp.

Reproduce: `./build/tests/test_cavity_eigenmodes`

Pass criterion: at least 6 of the first 8 analytical modes matched within
10%. This is a loose tolerance; lowest-order edge elements on the coarse
test mesh dominate the error.

## Lossy Dielectric Attenuation

Waveguide section filled with lossy dielectric, |S21| compared against the
analytical attenuation exp(-αL).

Reproduce: `PYTHONPATH=build/python:python python3 -m pytest python/tests/test_validation_suite.py -k lossy -v`

Pass criteria: measured attenuation within 35% of analytical, and lossy
|S21| strictly below lossless |S21|. Tolerance is wide because the lumped
port ABC interacts with the complex propagation constant.

## Dielectric Slab / Fresnel (analytical reference only)

`tests/test_dielectric_slab_fresnel.cpp` checks mesh setup (ports, material
regions) and verifies the analytical Fabry-Pérot reference satisfies
R + T = 1. **It does not compare an FEM solve against Fresnel** — that
requires Floquet ports, which are not implemented.
`examples/validation_fresnel.py` plots the analytical curves.

## What is NOT validated

- Radiation patterns against an analytical antenna (dipole/aperture); NTF
  tests are point-source sanity checks only.
- TEM / coax S-parameters (scalar port eigensolver cannot produce the TEM
  mode; the coax test is setup-sanity only).
- Periodic/unit-cell reflection against literature; oblique incidence.
- Dispersive materials inside an FEM solve (the material models themselves
  are unit-tested at the constitutive level).
- Lossy conductors, surface impedance, roughness.
- Behavior near or below waveguide cutoff.
- Any frequency outside roughly 5-26 GHz.

No comparison against commercial solvers (HFSS/CST) or measurement exists.

## Test Files

| Test | File | What it asserts |
|------|------|-----------------|
| WR-90 benchmark | `tests/benchmark_wr90.cpp` | S-params vs analytical, 7-12 GHz |
| WR-42 benchmark | `tests/test_wr42_waveguide.cpp` | S-params vs analytical, 20-26 GHz |
| ABC scaling | `tests/test_abc_scaling_sweep.cpp` | Optimal α in [0.4, 0.6] |
| Cavity modes | `tests/test_cavity_eigenmodes.cpp` | Eigenfrequencies vs f_mnp (10%) |
| Eigenmode S-params | `tests/test_eigenmode_sparams.cpp` | Transmission via eigenvector ports |
| Port eigenvector | `tests/test_eigenvector_waveport.cpp` | Diagnostic (no assertions) |
| Python suite | `python/tests/test_validation_suite.py` | Passivity, lossy attenuation |

## References

1. Pozar, D.M., *Microwave Engineering*, 4th Ed., Wiley, 2011.
2. Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd Ed., Wiley-IEEE, 2014.

## Continuous Validation

CI (`.github/workflows/ci.yml`) builds and runs the full registered CTest
suite (including the benchmarks above, label `benchmark`) plus
`pytest python/tests/` on every push and pull request.
