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
Port formulation: assembled surface-mass-matrix ABC with 2D discrete port
modes (July 2026); see "Port Formulation" below.

Reproduce: `./build/tests/benchmark_wr90`

| Freq (GHz) | \|S11\| | \|S21\| | Phase(S21) | Expected | Phase error |
|-----------|-------|-------|------------|----------|-------------|
| 7.0 | 0.019 | 0.9996 | -150.6° | -147.1° | 3.4° |
| 8.0 | 0.003 | 0.9997 | 80.8° | 84.8° | 4.0° |
| 9.0 | 0.010 | 0.9994 | -15.6° | -10.1° | 5.4° |
| 10.0 | 0.013 | 0.9992 | -100.7° | -93.3° | 7.4° |
| 11.0 | 0.041 | 0.9982 | 179.9° | -170.3° | 9.8°* |
| 12.0 | 0.053 | 0.9973 | 103.9° | 116.6° | 12.7° |

*Phase wrapping near ±180°.

Test pass criteria (deliberately loose for CI stability): |S11| < 0.15,
|S21| > 0.90, phase error < 15°, passivity ≤ 1.05. Honest summary at this
mesh density: **|S21| within 0.3%, |S11| below 0.06, phase error growing
to ~13° at the band edge** (phase error is discretization dispersion and
shrinks ~O(h²) with refinement — see Mesh Convergence).

A WR-42 benchmark (20-26 GHz) runs the same checks: `./build/tests/test_wr42_waveguide`.

## Port Formulation and the ABC Scale Factor

The wave-port boundary is a first-order modal ABC assembled as
`A += jβ·M_s`, where M_s is the port surface mass matrix
∫(n̂×N_i)·(n̂×N_j)dS, with the port mode taken from a 2D generalized
eigensolve (K_s, M_s) on the port face and β computed from the discrete
cutoff. The theoretical scale of the operator is exactly 1.0, and the
sweep confirms it — the former empirical `port_abc_scale = 0.5` belonged
to a diagonal-lumped approximation and is gone.

Reproduce: `./build/tests/test_abc_scaling_sweep` (writes `abc_scaling_sweep.csv`)

| α | \|S11\| | \|S21\| |
|-----|-------|-------|
| 0.80 | 0.207 | 0.978 |
| 0.90 | 0.094 | 0.995 |
| **1.00** | **0.013** | **0.999** |
| 1.20 | 0.187 | 0.982 |
| 1.50 | 0.387 | 0.922 |

Measured optimum: α = 1.00 (test asserts α_opt ∈ [0.9, 1.1] and
|S11| < 0.05 at α = 1.0). Remaining limitations: single-mode first-order
ABC (higher-order/evanescent content at the port is absorbed with the
dominant-mode impedance, not mode-matched), TE modes of hollow guides only.

## Mesh Convergence (WR-90 at 10 GHz)

Reproduce: `PYTHONPATH=build/python:python python3 scripts/run_convergence_study.py`

Measured (macOS ARM64, July 2026, assembled-M_s ports):

| Elements/λ | Tets | \|S21\| | \|S21\| error | \|S11\| | Phase error | Runtime |
|-----------|------|-------|-------------|-------|-----------|---------|
| 5 | 452 | 0.9868 | 1.32% | 0.150 | 19.1° | <0.1 s |
| 8 | 1,323 | 0.9978 | 0.22% | 0.037 | 8.3° | <0.1 s |
| 10 | 2,462 | 0.9995 | 0.05% | 0.007 | 4.1° | 0.1 s |
| 14 | 5,740 | 0.9996 | 0.04% | 0.006 | 2.2° | 0.6 s |
| 18 | 12,233 | 0.9999 | 0.01% | 0.0005 | 0.9° | 11.6 s |

All quantities now converge monotonically with refinement: phase error at
the expected ~O(h²) dispersion rate (19.1° → 0.9°), |S11| from 0.15 to
5×10⁻⁴, and |S21| error from 1.3% to 0.01% (fitted order ≈ 3.5 over this
range). Before the assembled-M_s port formulation, S-parameter magnitudes
plateaued at a 1-2% floor set by the diagonal-lumped ABC regardless of
mesh density; that plateau is eliminated.

Runtime is dominated by the sparse solve (the former dense 3D port
eigenvector computation — 1014 s at 12k tets — is replaced by a 2D dense
eigensolve on the port face, <0.1 s).

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
| ABC scaling | `tests/test_abc_scaling_sweep.cpp` | Optimal α in [0.9, 1.1]; |S11| < 0.05 at α=1 |
| Cavity modes | `tests/test_cavity_eigenmodes.cpp` | Eigenfrequencies vs f_mnp (10%) |
| Eigenmode S-params | `tests/test_eigenmode_sparams.cpp` | Transmission via eigenvector ports |
| Triangle edge mass | `tests/test_triangle_mass_matrix.cpp` | Closed form vs quadrature (1e-12) |
| Python suite | `python/tests/test_validation_suite.py` | Passivity, lossy attenuation |

## References

1. Pozar, D.M., *Microwave Engineering*, 4th Ed., Wiley, 2011.
2. Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd Ed., Wiley-IEEE, 2014.

## Continuous Validation

CI (`.github/workflows/ci.yml`) builds and runs the full registered CTest
suite (including the benchmarks above, label `benchmark`) plus
`pytest python/tests/` on every push and pull request.
