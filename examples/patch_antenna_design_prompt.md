# Single-shot prompt: Rectangular Microstrip Patch (2.45 GHz ISM)

**Goal:** Design and simulate a 2.45 GHz rectangular microstrip patch antenna on FR-4 with an inset 50 Ω microstrip feed. Validate S₁₁, input impedance, bandwidth, and far-field gain/pattern. Provide a lightweight param sweep + an optional auto-tune.

## Specs & targets

* Center frequency: **2.45 GHz**
* Return loss: **S₁₁ ≤ −10 dB** at 2.45 GHz; bandwidth goal ≥ **70–100 MHz**
* Realized broadside gain: **≥ 6 dBi** (FR-4 is lossy; accept ≥ 5.5 dBi)
* Efficiency: **≥ 50 %** (stretch ≥ 65 %)
* Polarization: linear (broadside)

## Stackup, materials, board size

* Substrate: **FR-4**, εr = **4.4**, tanδ = **0.02**, thickness **h = 1.6 mm**
* Copper: thickness **35 µm**, σ = **5.8e7 S/m** (use finite conductivity, not PEC)
* Board (finite ground) size: **80 mm (x) × 70 mm (y) × 1.6 mm (z)**
  Ground plane: full board underside (y=0 face). Patch sits on top (y=h).

## Patch geometry (initial closed-form seed at 2.45 GHz on FR-4, h=1.6 mm)

* Patch **width W ≈ 37.3 mm**, **length L ≈ 28.8 mm** (length is along the feed/inset direction)
* Place patch centered in x, with its radiating edge parallel to x.
* **Inset microstrip feed** centered along patch width (x-centered).
  Feed line width for 50 Ω on FR-4, h=1.6 mm: **W₅₀ ≈ 3.08 mm**.
* Cut an inset slot of width **W₅₀** into the patch from one radiating edge; **initial inset depth y₀ = 7.0 mm**.

## Feed & ports

* Implement a **microstrip feed line** (W₅₀ ≈ 3.08 mm) from the board edge to the inset gap; join to the patch at the end of the inset.
* Port: **50 Ω wave port** or **lumped port** between the feed line’s signal trace and ground reference.

  * If a lumped port: place across a short gap in the feed line near the board edge; define the port reference to ground.

## Surroundings & boundary condition

* Enclose the PCB in an **air box** that extends ≥ **λ₀/4 ≈ 30.6 mm** beyond all board edges and above the patch.
* Outer boundary: **PML** preferred; otherwise **radiation (absorbing) boundary**.

## Meshing & solve setup

* Curvilinear/tetrahedral FEM with:

  * Max element size in air: **λ₀/10** at 2.45 GHz
  * In dielectric: **λg/10** (auto-compute from εeff), and **edge refinement** on patch perimeters, feed, and inset (min size \~ **0.2 mm**)
  * **Adaptive meshing**: 3–5 passes or until |ΔS₁₁|<0.02 at 2.45 GHz
* Frequency sweep: **2.20–2.70 GHz**, coarse first (e.g., 10–20 points), then **adaptive/narrowband refinement** around the min |S₁₁|.

## Outputs (required)

1. **S-parameters:** Plot S₁₁(f); report f₍min|S11|₎, |S₁₁| at 2.45 GHz, and −10 dB bandwidth. Export **.s1p**.
2. **Input impedance:** R+jX vs f; value at 2.45 GHz.
3. **Fields:** |E| on a cut plane and surface currents at 2.45 GHz.
4. **Far-field:** 3D realized gain pattern at 2.45 GHz; 2D cuts (φ=0°, 90°). Report **peak realized gain**, **front-to-back**, **radiation/total efficiency**.
5. **Mesh report:** element count, min/avg edge length, adaptive convergence table.

## Quick param sweeps (exercise your param engine)

* **Inset depth y₀:** sweep **4 to 11 mm** in 1 mm steps; pick best S₁₁ at 2.45 GHz.
* **Patch length L:** local fine sweep **(L ± 0.6 mm)** in 0.2 mm steps (controls resonance).
* (Optional) **Loss model toggle:** tanδ ∈ {0.02, 0.015, 0.01} to show efficiency/gain sensitivity.

## One-click auto-tune (optional)

* Objective: minimize |S₁₁(2.45 GHz)| with variables {L, y₀} bounded as above; keep W fixed.
* Constraint: maintain −10 dB bandwidth ≥ 70 MHz if possible.
* Return tuned dims and re-run final solves + far-field.

## Acceptance criteria (pass/fail bands)

* **Match:** |S₁₁(2.45 GHz)| ≤ −10 dB (good: ≤ −15 dB)
* **Bandwidth (−10 dB):** ≥ 70 MHz (good: ≥ 90 MHz)
* **Realized gain (broadside):** ≥ 5.5 dBi (good: ≥ 6.5 dBi)
* **Efficiency:** ≥ 50 %
* **Pattern:** single broadside main lobe; cross-pol ≤ −15 dB at boresight; F/B ≥ 15 dB

## Exports

* Touchstone: **patch\_2p45\_FR4.s1p**
* CSV: S₁₁ vs f, Zin vs f
* PNG/SVG: S₁₁ plot, 2D cuts, surface currents
* VTK/OBJ: far-field mesh (if supported)
* JSON/YAML: final tuned geometry with units

---

## (Optional) Second substrate variant for regression

Re-run the exact workflow on **RO4350B** (εr = 3.66, tanδ = 0.0037, h = 1.524 mm).
Seed dimensions at 2.45 GHz: **W ≈ 40.1 mm**, **L ≈ 31.6 mm**, **50 Ω line width ≈ 3.36 mm**.
Use identical sweeps/metrics to compare bandwidth, gain, and efficiency vs FR-4.

---

### Notes for your simulator QA

* This case stresses: layered dielectrics, finite-conductivity metals, thin features, inset slot field singularities, PML stability, adaptive refinement on resonant geometry, and far-field integration.
* Keep geometry parametric so your UI can regenerate CAD from the final JSON/YAML.
* For reproducibility, log solver tolerances, mesh seeds, and adaptive stopping criteria with the results.

If you want, I can turn this into a ready-to-parse YAML (schema of your choosing) or a CLI script template that sets every knob exactly.
