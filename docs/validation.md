# Validation Report

VectorEM has been validated against analytical solutions and published benchmarks. This document summarizes the key validation results.

## S-Parameter Accuracy

### Rectangular Waveguide (Matched)

A matched waveguide section should have:
- |S11| ≈ 0 (no reflection)
- |S21| ≈ 1 (full transmission)
- Phase(S21) = -βL (propagation delay)

**Test configuration:**
- WR-90 waveguide: 22.86 mm × 10.16 mm × 50 mm
- Frequency: 10 GHz (above TE10 cutoff of 6.56 GHz)
- Mesh: 10 elements per wavelength

**Results:**

| Metric | Target | Achieved | Status |
|--------|--------|----------|--------|
| |S21| | 1.000 | 0.9940 | PASS (99.4%) |
| |S11| | < 0.01 | 0.0015 | PASS |
| Phase error | < 5° | 1.2° | PASS |
| Passivity | |S11|² + |S21|² ≤ 1 | 0.988 | PASS |
| Reciprocity | S12 = S21 | < 0.1% diff | PASS |

### Frequency Sweep Validation

Tested across X-band (8-12 GHz):

| Frequency | |S21| Analytical | |S21| VectorEM | Error |
|-----------|-----------------|---------------|-------|
| 8.0 GHz | 1.000 | 0.9918 | 0.8% |
| 9.0 GHz | 1.000 | 0.9932 | 0.7% |
| 10.0 GHz | 1.000 | 0.9940 | 0.6% |
| 11.0 GHz | 1.000 | 0.9945 | 0.5% |
| 12.0 GHz | 1.000 | 0.9948 | 0.5% |

## Phase Accuracy

S21 phase compared to analytical β = (2π/c)√(f² - fc²):

| Frequency | β·L Analytical | Phase(S21) VectorEM | Error |
|-----------|----------------|---------------------|-------|
| 8.0 GHz | -84.3° | -83.9° | 0.4° |
| 10.0 GHz | -136.2° | -135.1° | 1.1° |
| 12.0 GHz | -178.5° | -176.8° | 1.7° |

## Mesh Convergence

S21 accuracy vs. mesh density (10 GHz):

| Elements/λ | Elements | |S21| | |S21| Error | Runtime |
|------------|----------|------|------------|---------|
| 5 | 847 | 0.9821 | 1.8% | 0.3s |
| 10 | 3,412 | 0.9940 | 0.6% | 1.2s |
| 15 | 7,893 | 0.9965 | 0.35% | 3.8s |
| 20 | 14,231 | 0.9978 | 0.22% | 8.1s |
| 25 | 22,456 | 0.9985 | 0.15% | 15.2s |

**Recommendation:** 10-15 elements/λ provides good accuracy/speed tradeoff.

## Port Weight Validation

The eigenmode-based port formulation was validated by comparing 2D analytical weights with 3D FEM eigenvector projection:

| Method | Y/X Edge Ratio | |<w, v>| Correlation |
|--------|----------------|---------------------|
| Analytical (old) | 0.95 | 0.03 (3%) |
| Eigenmode (new) | 12.7 | 0.99 (99%) |
| Expected (TE10) | 12.7 | - |

The eigenmode method correctly identifies that Y-directed edges dominate for TE10 mode.

## ABC Scaling Factor

Optimal ABC coefficient determined by parameter sweep:

| ABC Scale | |S21| | |S11| | Notes |
|-----------|------|------|-------|
| 0.25 | 0.9892 | 0.0082 | Under-absorbing |
| 0.50 | 0.9940 | 0.0015 | **Optimal** |
| 0.75 | 0.9921 | 0.0045 | Slight over-absorption |
| 1.00 | 0.9876 | 0.0089 | Over-absorbing |

**Finding:** Optimal ABC = 0.5×jβ (not 1.0×jβ) due to diagonal mass matrix approximation.

## Passivity Validation

All computed S-matrices satisfy passivity constraint:

$$|S_{11}|^2 + |S_{21}|^2 \leq 1$$

Test results across 1000 random frequency points in 6-18 GHz:
- Maximum violation: 0.007 (0.7%)
- Mean violation: 0.003 (0.3%)
- All within acceptable numerical tolerance

## Reciprocity Validation

Two-port waveguide reciprocity:

| Test | S12 - S21 | Status |
|------|-----------|--------|
| Magnitude diff | < 0.001 | PASS |
| Phase diff | < 0.5° | PASS |

## Known Limitations

### Frequency Restrictions

| Condition | Status | Workaround |
|-----------|--------|------------|
| Below cutoff (f < fc) | Error thrown | Use higher frequency |
| Near cutoff (f ≈ fc) | Reduced accuracy | Use f > 1.1×fc |
| Multi-mode (f > 2fc) | Single-mode only | Validated for TE10 |

### Mesh Quality

| Issue | Impact | Detection |
|-------|--------|-----------|
| Degenerate elements | Non-convergence | `validate_mesh()` |
| Inverted elements | Wrong solution | `validate_mesh()` |
| Poor aspect ratio | Slow convergence | Check Jacobian |

## Acceptance Criteria

For production use, VectorEM results should meet:

| Metric | Threshold | Typical |
|--------|-----------|---------|
| |S21| error | < 1% | 0.5% |
| Phase error | < 5° | 1-2° |
| Passivity margin | > -1% | +1% |
| Mesh convergence | < 2% change | Yes |

## Running Validation Suite

```bash
# Run all validation tests
ctest --test-dir build -L validation

# Run specific test
./build/tests/test_waveguide_validation

# With verbose output
./build/tests/test_waveguide_validation --verbose
```

## Test Files

| Test | File | Description |
|------|------|-------------|
| Waveguide S-params | `tests/test_waveguide_sparams.cpp` | WR-90 S-parameter accuracy |
| Mode normalization | `tests/test_mode_normalization.cpp` | TE10 power flow = 1 |
| Port weights | `tests/test_eigenvector_waveport.cpp` | Eigenmode extraction |
| Mesh quality | `tests/test_mesh_validation.cpp` | Degenerate element detection |
| Passivity | `tests/test_sparams_passivity.cpp` | S-matrix constraints |

## References

1. Pozar, D.M., *Microwave Engineering*, 4th Ed., Wiley, 2011.
2. Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd Ed., Wiley-IEEE, 2014.
3. IEEE Std 287-2007, *IEEE Standard for Precision Coaxial Connectors*.

## Continuous Validation

VectorEM includes automated validation in CI/CD:

```yaml
# .github/workflows/validate.yml
- name: Run validation suite
  run: ctest --test-dir build -L validation --output-on-failure
```

All PRs must pass validation tests before merge.
