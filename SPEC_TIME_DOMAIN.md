/Users/johnhodge/Documents/code/VectorEM/system_design_spec.mdBelow is a proposed **research-grade specification** for a new **time-domain solver package** that complements your current EdgeFEM stack.

I’m grounding this in what EdgeFEM already is today: a **frequency-domain** 3D FEM solver for RF/mmWave using **Nédélec edge elements**, ([GitHub][1])ML/ABC boundaries, dispersive materials, Gmsh-based meshing, pybind11-backed Python bindings, and phased-array / metasurface workflows**. The current repo and PyPI package also state that **time-domain solving is not implemented in v1.0**. ([PyPI][2])

## Proposed product

**Working name:** `edgefem.td`
**Positioning:** a time-domain backend inside the EdgeFEM ecosystem, not a separate unrelated codebase.

That choice is deliberate. EdgeFEM already has the right outer architecture for this: Python SDK + C++ core, design classes, plotting/export utilities, phased-array examples, and an ecosystem link to higher-level phased-array modeling. A time-domain package should reuse that interface style instead of creating a parallel user experience. ([GitHub][1])

---

# 1. Purpose

Build a **broadband transient electromagnetic solver** for:

* **finite phased arrays**
* **finite metasurfaces**
* **wideband antennas and apertures**
* **transient near fields and far fields**
* **port-to-port impulse responses**
* **wideband active impedance / embedded element pattern analysis**

The package should target **PhD-level researchers and advanced RF/EM engineers**, not casual users.

## Primary value

Where EdgeFEM’s current frequency-domain engine is strongest for narrowband or per-frequency full-wave analysis, this package should be strongest when the user wants:

* one simulation for **many frequencies**
* **pulse / UWB** behavior
* **time delay steering**
* finite-aperture edge effects
* dispersive wideband substrates and radomes
* port/coupling behavior over a wide band in one run

---

# 2. Design philosophy

## Guiding principle

Do **not** try to clone CST or XFDTD in v1.

Instead, build the smallest serious solver that is:

* scientifically credible
* conformal to complex RF geometry
* usable from Python
* extensible for papers and prototypes
* tightly aligned with existing EdgeFEM abstractions

## Non-goals for v1

Do not include these in the first serious release:

* GUI-first workflow
* MPI distributed solve
* GPU backend
* hp-adaptivity
* full SPICE/circuit co-simulation
* topology optimization / adjoint inverse design
* nonlinear active device models
* full general CAD healing stack

Those are valuable, but they will slow the core solver down too early.

---

# 3. Numerical method recommendation

## Recommended core method

Use an **unstructured edge-element time-domain solver**, specifically:

* **first-order Maxwell equations**
* **curl-conforming edge basis**
* **DGTD-style explicit formulation** on tetrahedral/prismatic meshes

In plainer terms: build a **discontinuous Galerkin time-domain (DGTD)** solver using **Nédélec-like edge-oriented field representation**.

## Why this is the right choice

This is the best fit for your ecosystem because:

1. **Geometry fidelity matters**
   Finite phased arrays and metasurfaces have traces, vias, slots, tapered feeds, curved radomes, finite ground planes, and small discontinuities. A conformal unstructured mesh is a much better fit than a pure Yee-grid worldview.

2. **It matches EdgeFEM conceptually**
   EdgeFEM already centers edge-based EM discretization, Gmsh meshing, and RF-oriented port workflows. This preserves intellectual and software continuity. ([GitHub][1])

3. **Time-domain should avoid a large global linear solve at every time step**
   For broadband transient work, an explicit local formulation is much more attractive than implicit TD-FEM.

4. **It is research-grade**
   DGTD is a credible, publishable, extensible foundation for arrays, surfaces, dispersive media, and high-order extensions.

## Why not a Yee FDTD core

A Yee/FDTD solver is simpler, but I would not make it the primary architecture here because:

* it is weaker for conformal finite-array/metasurface geometry
* thin curved metal/dielectric features become a subcell-method problem
* meshing of vias, probes, edge chamfers, and finite radomes is less natural
* it drifts away from the design language of EdgeFEM

FDTD could be a future lightweight backend, but not the flagship one.

## Why not implicit TD-FEM

Implicit TD-FEM is attractive mathematically, but for this package it is the wrong first choice because:

* each time step requires a large solve
* that is painful for long transients
* broadband arrays need many steps
* implementation complexity rises quickly

---

# 4. Governing equations

Use first-order Maxwell in differential form:

[
\frac{\partial \mathbf{B}}{\partial t} = - \nabla \times \mathbf{E}
]

[
\frac{\partial \mathbf{D}}{\partial t} = \nabla \times \mathbf{H} - \mathbf{J}
]

with constitutive relations

[
\mathbf{D} = \varepsilon_0 \varepsilon_r \mathbf{E} + \sum_k \mathbf{P}_k
]

[
\mathbf{B} = \mu_0 \mu_r \mathbf{H}
]

where (\mathbf{P}_k) are dispersive polarization states.

## Material update model

Use **ADEs** as the primary dispersion framework:

* Debye
* Lorentz
* Drude
* Drude-Lorentz

That is consistent with the material scope already exposed by EdgeFEM today. ([PyPI][2])

---

# 5. Scope of supported problems

## v1.0 must support

### Finite phased arrays

* linear, planar, and sparse arrays
* element-level port excitation
* simultaneous multiport pulse excitation
* phase steering and true-time-delay steering
* embedded element patterns
* active impedance vs scan angle
* finite-array edge effects
* mutual coupling over frequency

### Finite metasurfaces

* finite apertures built from repeated or nonuniform cells
* reflective and transmissive sheets
* beam deflection / focusing profiles
* dispersive substrate-backed patterned layers
* finite-size diffraction and truncation effects

### Antennas and wideband radiators

* probe-fed and waveguide-fed structures
* patch, slot, horn, loaded aperture, conformal radiator classes
* transient radiation and time-domain near-field probes

---

# 6. Package architecture

## Recommendation

Keep the repo unified for now, with a new namespace:

* `edgefem.fd` for the existing frequency-domain solver
* `edgefem.td` for the new time-domain solver

Do **not** split into a separate repo yet.

That avoids premature fragmentation and lets you reuse:

* Gmsh geometry / tagging workflows
* material definitions
* port definitions
* design classes
* plotting/export layer
* docs and CI conventions

EdgeFEM already has a clear Python/C++ split, a C++ core, and a research-oriented roadmap. This is the right place to extend, not fork. ([GitHub][1])

## Core modules

### Shared/core

* mesh IO and tags
* units/constants
* material library
* geometry/design classes
* boundary-condition schemas
* common export formats

### New TD modules

* `td/operators`
* `td/materials`
* `td/pml`
* `td/ports`
* `td/sources`
* `td/timesteppers`
* `td/probes`
* `td/ntff`
* `td/post`
* `td/array`
* `td/checkpoint`

## Language/runtime

* **C++20 core**
* **pybind11 bindings**
* Python 3.9+ user API

That matches the current package style. ([PyPI][2])

---

# 7. Discretization details

## Mesh types

Support:

* tetrahedral volume mesh
* prism layers where useful for stacked dielectrics and ground-backed structures
* triangular surface meshes for impedance or GSTC sheets

## Basis functions

### v1

* first-order curl-conforming edge basis in volume
* consistent tangential-trace basis on faces
* first-order surface basis for sheet operators

### v1.1+

* higher-order basis
* curved elements if practical

## Numerical flux

Use standard DG numerical fluxes with selectable modes:

* centered flux
* upwind flux
* low-reflection tuned flux for radiating problems

Default should be **upwind**.

---

# 8. Time integration

## Default

Use **low-storage explicit Runge-Kutta (LSERK4)** as the default production integrator.

## Optional fast path

Provide **leapfrog / staggered explicit stepping** for simple nondispersive problems.

## Why

LSERK4 is the safer default because it handles:

* dispersive ADE coupling
* source injection cleanly
* more general operator splitting
* future local-time-stepping compatibility

## CFL handling

The engine must compute and expose:

* global stable time step
* element-wise CFL bottlenecks
* per-region min time step estimates

Users need to know when tiny mesh features are destroying runtime.

---

# 9. Boundary conditions

## v1 required

* PEC
* PMC
* symmetry planes
* first-order absorbing boundary
* **CPML / UPML**
* radiation box surfaces for NTFF

## v1.1 required

* Bloch-periodic boundaries
* phase-shifted periodic boundaries
* periodic benchmarking mode for unit-cell cross-checks

For finite metasurfaces, periodic boundaries are not the main production feature, but they are very useful for validating finite-aperture behavior against unit-cell baselines.

---

# 10. Source and port models

## Port models

### v1

* lumped voltage/current port
* coaxial feed port
* waveguide modal port
* differential port pairs
* multiport simultaneous excitation

### v1.1

* broadband modal de-embedding
* port calibration/reference plane shift
* active reflection extraction under steering law

## Source models

* Gaussian pulse
* differentiated Gaussian
* modulated sinusoid
* Ricker wavelet
* user-defined waveform
* array excitation table: amplitude, phase, delay per element

## Future

* total-field / scattered-field plane-wave injection
* Huygens surface source
* impressed equivalent currents from aperture synthesis

---

# 11. Material model

## v1 required

* isotropic lossless dielectric
* isotropic lossy dielectric
* conductivity
* Debye/Lorentz/Drude/Drude-Lorentz via ADE
* PEC as idealized metal

## Strongly recommended for metasurfaces

Add **surface models** early:

* surface impedance
* surface admittance
* resistive sheet
* thin conductive coating model

This is crucial because many metasurface problems become unnecessarily expensive if every ultrathin layer must be volumetrically meshed.

## v1.1

* diagonal anisotropy
* magnetic loss
* tensor surface impedance

## v2

* full tensor anisotropy
* bianisotropic sheet/GSTC operators

---

# 12. Finite-array and metasurface-first workflows

This package should not look like a generic EM code with arrays as an afterthought.

It needs dedicated workflows.

## Array workflow objects

* `ArrayLattice`
* `ElementInstance`
* `ExcitationLaw`
* `ScanSchedule`
* `EmbeddedPatternSet`
* `ActiveImpedanceReport`

## Metasurface workflow objects

* `SurfaceTile`
* `ApertureLayout`
* `PhaseProfile`
* `SheetModel`
* `TransmissionReflectionReport`

## Required outputs

* port voltages/currents vs time
* broadband S-parameters from FFT
* embedded element patterns
* active S/Z under steering
* realized gain vs frequency and scan
* beam pointing error
* beam squint
* sidelobe level
* aperture phase/amplitude maps
* finite-aperture truncation loss
* near-field movies
* probe histories
* wideband NTFF patterns

---

# 13. Near-to-far-field and postprocessing

## Near-to-far transform

Use a closed-surface NTFF implementation around the radiator / aperture region.

## Wideband far field

One transient run should support postprocessed far fields at many frequencies.

## Postprocessing API

The solver should expose:

* FFT windowing choices
* gating
* time-domain waveform visualization
* frequency slicing
* polarization decomposition
* array manifold export
* pattern cuts and 3D patterns

The current EdgeFEM package already has plotting and export conventions for patterns, coupling, Touchstone, VTK, and CSV. Reuse those formats wherever possible. ([GitHub][1])

---

# 14. Proposed Python API

```python
from edgefem.td import (
    TransientSimulation,
    ArrayDesign,
    MetasurfaceDesign,
    GaussianPulse,
    LSERK4,
    CPML,
)

arr = ArrayDesign.rectangular_patch_array(
    nx=8, ny=8,
    dx=0.55, dy=0.55,              # wavelengths or SI by config
    substrate_eps_r=2.2,
    substrate_h=1.6e-3,
    ground_plane_margin=10e-3,
)

arr.assign_feed_ports(kind="coax")
arr.generate_mesh(density=18)

sim = TransientSimulation(
    design=arr,
    solver=LSERK4(),
    boundary=CPML(thickness_cells=10),
)

sim.set_excitation(
    ports="all",
    waveform=GaussianPulse(f0=10e9, fbw=0.6),
    steering={"theta_deg": 25, "phi_deg": 0, "mode": "time_delay"},
)

sim.add_near_to_far_box()
sim.add_probe(point=(0.0, 0.0, 20e-3), field="E")
sim.run(t_end=12e-9)

S = sim.sparams(freqs=[8e9, 9e9, 10e9, 11e9, 12e9])
ff = sim.farfield(freqs=[9e9, 10e9, 11e9])
Zact = sim.active_impedance(scan_theta=[0, 10, 20, 30], phi_deg=0)
```

That API should feel like EdgeFEM, not like a low-level PDE lab notebook.

---

# 15. Performance model

## v1 target platform

Optimize for:

* high-end workstation
* shared-memory CPU
* 64–256 GB RAM
* OpenMP parallel loops
* restart/checkpoint support

## v1 optimization priorities

1. contiguous memory layout
2. element-block batching
3. sparse/local operator reuse
4. low-overhead probe recording
5. efficient FFT/postprocessing
6. checkpoint/restart for long runs

## Not in v1

* GPU
* MPI
* task-based distributed decomposition

Those belong later, consistent with the broader staged roadmap style already present in EdgeFEM. ([GitHub][1])

---

# 16. Verification and validation plan

This part matters a lot. A PhD-level solver lives or dies by V&V.

## Tier 1: mathematical / code verification

* manufactured solution test
* curl-operator consistency
* energy conservation in lossless closed problems
* convergence with mesh refinement
* convergence with time-step refinement
* CPML reflection benchmark
* reciprocity checks where applicable

## Tier 2: canonical EM benchmarks

* 1D pulse propagation in homogeneous medium
* rectangular cavity transient response
* rectangular waveguide pulse transmission vs analytic modal solution
* coax or lumped-fed strip/patch benchmark
* dipole / patch far-field benchmark

## Tier 3: EdgeFEM cross-validation

Use the existing frequency-domain engine as a direct reference for:

* waveguide S-parameters
* patch resonance / input impedance
* unit-cell reflection/transmission
* phased-array active impedance scan cuts

That is especially valuable because EdgeFEM already exposes waveguide, patch, unit-cell, and phased-array workflows. ([GitHub][1])

## Tier 4: application validation

* finite 4x4 and 8x8 array beam steering
* scan loss vs angle
* embedded pattern extraction
* finite metasurface anomalous reflection/transmission
* finite aperture edge diffraction and truncation

## Suggested acceptance criteria

* amplitude error within engineering tolerance across band
* phase error small enough for array steering use
* passivity preserved in passive structures
* no late-time PML instability
* NTFF main beam agreement within a defined dB / angle tolerance

---

# 17. Documentation deliverables

## Must-have docs

* theory manual
* developer architecture note
* user quickstart
* phased-array tutorial
* metasurface tutorial
* dispersion/material tutorial
* validation report
* API reference
* performance guide
* limitations page

EdgeFEM already has a docs/tutorial/validation pattern. Mirror that exactly. ([GitHub][1])

---

# 18. v1 roadmap

## v0.1 research prototype

* tetra mesh import
* basic DGTD operators
* PEC + simple absorber
* lumped source
* probe recording
* 1D/2D/3D smoke tests

## v0.5 engineering alpha

* CPML
* dispersive ADE
* modal/lumped/coax ports
* FFT-based S-parameters
* NTFF
* patch and waveguide validation

## v1.0 first serious release

* finite phased-array workflow
* time-delay steering
* active impedance extraction
* embedded patterns
* finite metasurface aperture workflow
* sheet impedance models
* checkpoint/restart
* publication-grade validation package

## v1.1

* Bloch periodic mode
* diagonal anisotropy
* higher-order basis
* better port calibration
* performance tuning

## v2.0

* local time stepping
* MPI
* GPU
* tensor/bianisotropic sheets
* optimization hooks

---

# 19. Recommendation on package boundaries

My recommendation is:

* **Keep this inside the EdgeFEM repo initially**
* expose it as **`edgefem.td`**
* share mesh/material/IO/design abstractions
* keep the solver kernels separate from the existing frequency-domain solver

That gives you reuse without forcing everything into one backend.

A separate package is only worth it later if:

* release cadence diverges sharply
* dependencies diverge sharply
* contributors split into separate teams
* code size becomes hard to maintain in one repo

Right now, splitting would be more organizational overhead than value.

---

# 20. Bottom-line recommendation

If I were writing the actual implementation plan, I would build:

* **an unstructured explicit DGTD solver**
* **inside the current EdgeFEM ecosystem**
* with **finite-array and metasurface workflows as first-class objects**
* and with **surface-sheet models added early** to avoid absurd thin-layer meshing costs

That is the best balance of:

* research credibility
* reuse of your current EdgeFEM architecture
* relevance to finite phased arrays and metasurfaces
* and realistic implementation scope

If you want, I can turn this into a proper **`SYSTEM_DESIGN_SPEC.md`** with sections like requirements, architecture diagrams, class layout, milestones, and acceptance tests.

[1]: https://github.com/jman4162/EdgeFEM "GitHub - jman4162/EdgeFEM: 3D finite-element electromagnetics (FEM) solver for RF/mmWave simulation. Computes S-parameters, radiation patterns, and fields using Nédélec edge elements. Supports wave ports, PML/ABC boundaries, periodic structures, and dispersive materials. · GitHub"
[2]: https://pypi.org/project/edgefem/1.0.0/ "edgefem · PyPI"
