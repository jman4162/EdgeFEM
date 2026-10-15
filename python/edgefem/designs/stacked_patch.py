"""
Stacked Patch Antenna Design Class

Multi-layer stacked patch antenna with a probe feed and a lumped port in
the probe gap at the ground plane.
Uses full-wave FEM simulation via lumped port excitation and
Stratton-Chu near-to-far field transformation for radiation patterns.
"""

import subprocess
import tempfile
import os
import types
import warnings
from typing import Optional, Tuple, List, Dict, Sequence, Union
import numpy as np

import pyedgefem as em


class _AperturePattern:
    """Lightweight far-field pattern compatible with FFPattern3D interface."""

    def __init__(self, theta_grid, phi_grid, E_theta, E_phi):
        self.theta_grid = theta_grid
        self.phi_grid = phi_grid
        self.E_theta = E_theta
        self.E_phi = E_phi

    def power_pattern(self):
        return np.abs(self.E_theta)**2 + np.abs(self.E_phi)**2

    def total_magnitude(self):
        return np.sqrt(np.abs(self.E_theta)**2 + np.abs(self.E_phi)**2)


class FrequencySweepResult(tuple):
    """(frequencies, Z_in, S11) from frequency_sweep(), plus the frequency
    whose field keep_fields="min_s11" kept.

    Unpacks as the three arrays, as before. min_s11_index (position in
    frequencies) and min_s11_frequency (Hz) are None unless keep_fields
    asked for "min_s11" and some frequency had a finite S11.
    """

    def __new__(cls, freqs, z_in, s11, min_s11_index=None):
        self = super().__new__(cls, (freqs, z_in, s11))
        self.min_s11_index = (None if min_s11_index is None
                              else int(min_s11_index))
        self.min_s11_frequency = (None if min_s11_index is None
                                  else float(freqs[min_s11_index]))
        return self

    def __getnewargs__(self):
        # pickle and copy rebuild through __new__
        return (*self, self.min_s11_index)


class StackedPatchDesign:
    """
    Multi-layer stacked patch antenna with a probe feed.

    The antenna consists of one or more stacked patches on dielectric
    substrates above a ground plane. A probe runs from the ground plane
    to the first patch; a lumped port bridges the gap between the probe
    and the ground plane. Full-wave FEM simulation provides
    S-parameters, input impedance, and radiation patterns.

    Geometry (bottom to top):
        - Ground plane (PEC) at z = 0, the bottom of the domain, with a
          square hole around the probe. The port strip spans the hole
          from the probe to the ground on the +x side; the rest of the
          hole is left without a boundary condition.
        - Probe: a square PEC post from z = 0 to the first patch, cut out
          of the mesh. Its side is probe_radius / 0.59, the thin-wire
          equivalent radius of a square (Balanis, Antenna Theory,
          equivalent-radius table), so the post has the inductance of a
          round probe of probe_radius. It is square so that the port
          strip touches the probe along a straight edge and stays a
          rectangle (the lumped port's strip width W = area / gap length
          is exact only for rectangles).
        - One box per substrate layer (tags 110, 111, ...), the full
          ground-plane area each.
        - Patch conductors (PEC) inside the patch footprints.
        - Air box with ABC on its top and sides, at least lambda/2 above
          the top metal at the lowest frequency given to generate_mesh().
        - Huygens surface for NTF, halfway up the air box.

    There is no cavity below the ground plane. A coax-fed model with the
    coax below the ground needs a wave port on the coax cross-section,
    not a lumped strip.
    TODO(coax feed): build it with a wave port; see the EdgeFEM
    investigation of 2026-10-12, sections A3 and E.2.

    Example:
        design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
            probe_radius=0.65e-3,
            probe_offset=(0.0, -6e-3),
        )
        design.generate_mesh(density=12, design_freq=2.4e9, f_min=2.2e9)
        S11, solution = design.simulate(2.4e9)
        pattern = design.radiation_pattern(2.4e9)
    """

    # Physical tag scheme
    _TAG_GROUND = 1
    _TAG_PROBE = 4
    _TAG_PORT = 5
    _TAG_PATCH_BASE = 10      # 10, 11, 12, ...
    _TAG_ABC = 50
    _TAG_HUYGENS = 60
    _TAG_SUBSTRATE_BASE = 110  # 110, 111, ...
    _TAG_AIR = 150

    # Thin-wire equivalent radius of a square cross-section, per unit side.
    _SQUARE_EQUIV_RADIUS = 0.59
    # Default hole/probe radius ratio: a 50 ohm air coax, ln(b/a) = 50/60.
    _DEFAULT_HOLE_RATIO = float(np.exp(50.0 / 60.0))
    # Default air box height above the top metal, in wavelengths at f_min.
    _AIR_BOX_WAVELENGTHS = 0.5
    # Default f_min as a fraction of design_freq: a sweep of +/-20% around
    # design_freq (the default band of PatchAntennaDesign.bandwidth) needs
    # no f_min and raises no air box warning.
    _F_MIN_DEFAULT_RATIO = 0.8
    # Mesh size field (see _set_mesh_size). Near the patches the size is
    # _PATCH_SIZE_X_DENSITY * (thinnest substrate thickness) / density:
    # h/3 at density 12. Chosen from the convergence study in
    # docs/validation/patch_2p45ghz.md (tools/patch_convergence.py): on the
    # 2.45 GHz FR-4 patch at density 12, 6.0 (h/2) gave f_res 2.360 GHz and
    # 4.0 gave 2.377 GHz, against 2.389 GHz from openEMS, at 3x the solve
    # time (97 s -> 303 s). The known-answer test in
    # python/tests/test_patch_physics.py depends on this value.
    _PATCH_SIZE_X_DENSITY = 4.0
    _PATCH_SIZE_SLOPE = 0.5
    _PROBE_SIZE_SLOPE = 0.5

    def __init__(
        self,
        patches: List[Dict],
        substrates: List[Dict],
        cavity_size: Optional[float] = None,
        cavity_depth: Optional[float] = None,
        probe_radius: float = 0.5e-3,
        probe_offset: Tuple[float, float] = (0.0, 0.0),
        z0: float = 50.0,
        ground_plane_size: Optional[float] = None,
        air_box_height: Optional[float] = None,
        *,
        probe_hole_radius: Optional[float] = None,
    ):
        """
        Initialize stacked patch design.

        Args:
            patches: List of patch dicts, each with keys:
                     'width' (m), 'length' (m), 'height' (m above ground)
            substrates: List of substrate dicts, each with keys:
                       'thickness' (m), 'eps_r', 'tan_d' (optional, default 0)
            cavity_size: Deprecated and ignored (the model has no cavity).
            cavity_depth: Deprecated and ignored (the model has no cavity).
            probe_radius: Radius of the probe feed (m). The probe is a
                square post with this equivalent radius.
            probe_offset: (x, y) offset of probe from center (m)
            z0: Reference impedance (ohms)
            ground_plane_size: Ground plane size. If None, 3x largest patch.
            air_box_height: Air box height above the top metal (m). If
                None, lambda/2 at the f_min given to generate_mesh().
            probe_hole_radius: Equivalent radius of the ground-plane hole
                around the probe (m). If None, 2.30 x probe_radius (the
                outer/inner ratio of a 50 ohm air coax).
        """
        if len(patches) == 0:
            raise ValueError("At least one patch is required")
        if len(substrates) == 0:
            raise ValueError("At least one substrate is required")
        if cavity_size is not None or cavity_depth is not None:
            warnings.warn(
                "StackedPatchDesign: cavity_size and cavity_depth are "
                "deprecated and ignored. The model has no cavity below the "
                "ground plane: the domain starts at z = 0 and the lumped "
                "port sits in the probe gap at the ground plane.",
                DeprecationWarning,
                stacklevel=2,
            )
        if probe_radius <= 0:
            raise ValueError("probe_radius must be positive")
        if probe_hole_radius is None:
            probe_hole_radius = self._DEFAULT_HOLE_RATIO * probe_radius
        if probe_hole_radius <= probe_radius:
            raise ValueError(
                f"probe_hole_radius ({probe_hole_radius}) must exceed "
                f"probe_radius ({probe_radius})"
            )

        self.patches = patches
        self.substrates = substrates
        self.cavity_size = cavity_size
        self.cavity_depth = cavity_depth
        self.probe_radius = probe_radius
        self.probe_hole_radius = probe_hole_radius
        self.probe_offset = probe_offset
        self.z0 = z0

        max_patch = max(max(p['width'], p['length']) for p in patches)
        self.ground_plane_size = ground_plane_size or 3.0 * max_patch
        self._air_box_height_arg = air_box_height
        # Set by generate_mesh(): the air box height used for the current
        # mesh, and the lowest frequency at which it is >= lambda/2.
        self._air_box_height = None
        self._air_box_f_min = None
        self._air_box_warned = False

        self._mesh = None
        self._mesh_path = None
        self._bc = None
        self._ports = None
        # Field solution per frequency. Written only by _store_solution().
        self._solutions = {}
        # Far-field cache: (freq, n_theta, n_phi) -> (solution, pattern),
        # least recently used first. An entry is used only while its
        # solution is still the one stored for freq.
        self._pattern_cache = {}

    @property
    def mesh(self):
        return self._mesh

    @property
    def _last_solutions(self):
        """Read-only view of the stored field solutions, keyed by frequency
        (Hz). Use _store_solution() to add or replace one."""
        return types.MappingProxyType(self._solutions)

    @property
    def air_box_height(self) -> Optional[float]:
        """Air box height above the top metal in the current mesh (m)."""
        return self._air_box_height

    def _probe_half_widths(self) -> Tuple[float, float]:
        """Half-sides (m) of the square probe and of the square hole."""
        k = 2.0 * self._SQUARE_EQUIV_RADIUS
        return self.probe_radius / k, self.probe_hole_radius / k

    def generate_mesh(
        self,
        density: float = 12.0,
        output_path: Optional[str] = None,
        design_freq: Optional[float] = None,
        f_min: Optional[float] = None,
    ) -> None:
        """
        Generate FEM mesh for the stacked patch antenna.

        Uses the Gmsh Python API with OpenCASCADE kernel for robust
        boolean operations on the complex antenna geometry.

        Args:
            density: Elements per wavelength
            output_path: Path for mesh file. If None, uses temp directory.
            design_freq: Design frequency for mesh sizing (Hz).
                        If None, estimates from first patch dimensions.
            f_min: Lowest frequency to be simulated on this mesh (Hz).
                Sets the default air box height to lambda/2 at f_min.
                If None, 0.8 * design_freq. simulate() and
                frequency_sweep() warn below the frequency at which the
                air box is lambda/2.
        """
        try:
            import gmsh
        except ImportError:
            raise RuntimeError(
                "Gmsh Python module required. Install with: pip install gmsh"
            )

        c0 = 299792458.0

        if design_freq is None:
            L = self.patches[0]['length']
            eps_r = self.substrates[0]['eps_r']
            design_freq = c0 / (2 * L * np.sqrt(eps_r))
        if f_min is None:
            f_min = self._F_MIN_DEFAULT_RATIO * design_freq

        wavelength = c0 / design_freq
        lc = wavelength / density

        if self._air_box_height_arg is not None:
            h_air = self._air_box_height_arg
        else:
            h_air = self._AIR_BOX_WAVELENGTHS * c0 / f_min
        self._air_box_height = h_air
        self._air_box_f_min = self._AIR_BOX_WAVELENGTHS * c0 / h_air
        self._air_box_warned = False

        if output_path is None:
            fd, msh_path = tempfile.mkstemp(suffix=".msh")
            os.close(fd)
        elif output_path.endswith(".msh"):
            msh_path = output_path
        else:
            msh_path = output_path + ".msh"

        self._build_gmsh_geometry(gmsh, lc, msh_path, density=density)

        self._mesh = em.load_gmsh(msh_path)
        self._mesh_path = msh_path
        self._bc = None
        self._ports = None
        self._solutions = {}
        self._clear_pattern_cache()

    # Patterns kept in the far-field cache. A 91 x 73 pattern holds about
    # 210 kB of arrays.
    _PATTERN_CACHE_MAX = 8

    def _solution_at(self, freq: float):
        """Field at freq: the stored one, else from simulate(freq)."""
        if freq not in self._solutions:
            self.simulate(freq)
        if freq not in self._solutions:
            raise RuntimeError(
                f"the solve at {freq / 1e9:.6g} GHz failed; no field")
        return self._solutions[freq]

    def _store_solution(self, freq: float, solution) -> None:
        """Record the field solution for freq and drop the cached patterns
        computed from the solution it replaces."""
        self._solutions[freq] = solution
        for key in [k for k in self._pattern_cache if k[0] == freq]:
            del self._pattern_cache[key]

    def _clear_pattern_cache(self) -> None:
        self._pattern_cache = {}

    def _check_air_box(self, freq: float) -> None:
        """Warn once per mesh when the air box is below lambda/2 at freq."""
        if self._air_box_f_min is None or self._air_box_warned:
            return
        if freq < self._air_box_f_min * (1.0 - 1e-9):
            c0 = 299792458.0
            warnings.warn(
                f"StackedPatchDesign: the air box is "
                f"{self._air_box_height * 1e3:.1f} mm above the top metal, "
                f"less than lambda/2 = {0.5 * c0 / freq * 1e3:.1f} mm at "
                f"{freq / 1e9:.4g} GHz; pass f_min={freq:.4g} to "
                f"generate_mesh() (or a larger air_box_height) to size it "
                f"for this frequency",
                stacklevel=3,
            )
            self._air_box_warned = True

    def _build_gmsh_geometry(self, gmsh, lc: float, msh_path: str,
                             density: float = 12.0):
        """Build antenna geometry using Gmsh Python API with OCC kernel."""
        gp = self.ground_plane_size
        h_air = self._air_box_height
        a, b = self._probe_half_widths()
        cx = gp / 2 + self.probe_offset[0]
        cy = gp / 2 + self.probe_offset[1]

        z_levels = [0.0]
        for s in self.substrates:
            z_levels.append(z_levels[-1] + s['thickness'])
        z_metal = max(z_levels[-1], max(p['height'] for p in self.patches))
        z_top = z_metal + h_air
        z_huygens = z_metal + h_air * 0.5  # Midpoint of air box
        z_probe = self.patches[0]['height']

        # The probe must end inside the first patch footprint, and the hole
        # must lie inside the ground plane.
        p0 = self.patches[0]
        if (abs(self.probe_offset[0]) + a > p0['width'] / 2
                or abs(self.probe_offset[1]) + a > p0['length'] / 2):
            raise ValueError(
                "probe (half-side {:.3g} m at offset {}) does not fit inside "
                "the first patch".format(a, self.probe_offset))
        if min(cx, cy) - b <= 0 or max(cx, cy) + b >= gp:
            raise ValueError(
                "probe hole (half-side {:.3g} m at offset {}) does not fit "
                "inside the ground plane".format(b, self.probe_offset))

        min_t = min(s['thickness'] for s in self.substrates)
        tol = 1e-3 * min(a, min_t)

        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Verbosity", 1)
            gmsh.model.add("stacked_patch")
            occ = gmsh.model.occ

            # Volumes: one box per substrate layer and the air box above,
            # with the probe post removed from all of them.
            boxes = []
            for i, s in enumerate(self.substrates):
                boxes.append((3, occ.addBox(0, 0, z_levels[i], gp, gp,
                                            s['thickness'])))
            boxes.append((3, occ.addBox(0, 0, z_levels[-1], gp, gp,
                                        z_top - z_levels[-1])))
            probe_box = occ.addBox(cx - a, cy - a, 0, 2 * a, 2 * a, z_probe)
            vols, _ = occ.cut(boxes, [(3, probe_box)])

            # z = 0: the hole around the probe (a square ring) and the port
            # strip in it, from the probe face at x = cx + a to the hole
            # edge at x = cx + b. E points along +x across the strip.
            hole = occ.addRectangle(cx - b, cy - b, 0, 2 * b, 2 * b)
            probe_foot = occ.addRectangle(cx - a, cy - a, 0, 2 * a, 2 * a)
            ring, _ = occ.cut([(2, hole)], [(2, probe_foot)])
            port = occ.addRectangle(cx + a, cy - a, 0, b - a, 2 * a)

            tools = list(ring) + [(2, port)]
            for patch in self.patches:
                tools.append((2, occ.addRectangle(
                    gp / 2 - patch['width'] / 2, gp / 2 - patch['length'] / 2,
                    patch['height'], patch['width'], patch['length'])))

            # Huygens surface for near-to-far field transformation, inset
            # from the domain boundary to avoid ABC-contaminated fields
            h_inset = gp * 0.1
            tools.append((2, occ.addRectangle(
                h_inset, h_inset, z_huygens,
                gp - 2 * h_inset, gp - 2 * h_inset)))

            occ.fragment(vols, tools)
            occ.synchronize()

            def inside(bb, x0, y0, x1, y1):
                return (bb[0] > x0 - tol and bb[1] > y0 - tol
                        and bb[3] < x1 + tol and bb[4] < y1 + tol)

            gnd_surfs, port_surfs, probe_surfs = [], [], []
            abc_surfs, huygens_surfs = [], []
            patch_surfs = {pi: [] for pi in range(len(self.patches))}
            for _, tag in gmsh.model.getEntities(2):
                bb = gmsh.model.getBoundingBox(2, tag)
                xe, ye, ze = bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2]
                zc = (bb[2] + bb[5]) / 2
                if ze < tol:  # horizontal
                    if abs(zc) < tol:
                        if inside(bb, cx + a, cy - a, cx + b, cy + a):
                            port_surfs.append(tag)
                        elif inside(bb, cx - b, cy - b, cx + b, cy + b):
                            pass  # rest of the probe hole: no BC
                        else:
                            gnd_surfs.append(tag)
                    elif abs(zc - z_top) < tol:
                        abc_surfs.append(tag)
                    elif abs(zc - z_huygens) < tol:
                        huygens_surfs.append(tag)
                    else:
                        # Patch PEC only inside each patch footprint; the
                        # rest of a substrate interface carries no BC.
                        for pi, patch in enumerate(self.patches):
                            if abs(zc - patch['height']) < tol and inside(
                                    bb,
                                    gp / 2 - patch['width'] / 2,
                                    gp / 2 - patch['length'] / 2,
                                    gp / 2 + patch['width'] / 2,
                                    gp / 2 + patch['length'] / 2):
                                patch_surfs[pi].append(tag)
                                break
                elif xe < tol or ye < tol:  # vertical
                    on_side = (
                        (xe < tol and (abs(bb[0]) < tol
                                       or abs(bb[0] - gp) < tol))
                        or (ye < tol and (abs(bb[1]) < tol
                                          or abs(bb[1] - gp) < tol)))
                    if on_side:
                        abc_surfs.append(tag)
                    elif inside(bb, cx - a, cy - a, cx + a, cy + a):
                        probe_surfs.append(tag)
                    else:
                        raise RuntimeError(
                            f"unclassified vertical surface {tag}, bbox {bb}")
                else:
                    raise RuntimeError(
                        f"unclassified surface {tag}, bbox {bb}")

            port_area = sum(occ.getMass(2, t) for t in port_surfs)
            if abs(port_area - (b - a) * 2 * a) > 1e-6 * (b - a) * 2 * a:
                raise RuntimeError(
                    f"port strip area {port_area:.6g} m^2, expected "
                    f"{(b - a) * 2 * a:.6g} m^2 (rectangle (b-a) x 2a)")

            groups = [
                (gnd_surfs, self._TAG_GROUND, "ground"),
                (port_surfs, self._TAG_PORT, "port"),
                (probe_surfs, self._TAG_PROBE, "probe"),
                (abc_surfs, self._TAG_ABC, "ABC"),
                (huygens_surfs, self._TAG_HUYGENS, "Huygens"),
            ]
            for pi, sfs in patch_surfs.items():
                groups.append((sfs, self._TAG_PATCH_BASE + pi, f"patch {pi}"))
            for sfs, ptag, name in groups:
                if not sfs:
                    raise RuntimeError(f"no surfaces found for the {name}")
                gmsh.model.addPhysicalGroup(2, sfs, ptag)

            # Volumes: substrate layer by the z-centroid, air above.
            substrate_vols = {si: [] for si in range(len(self.substrates))}
            air_vols = []
            for _, vtag in gmsh.model.getEntities(3):
                zc = occ.getCenterOfMass(3, vtag)[2]
                for si in range(len(self.substrates)):
                    if z_levels[si] < zc < z_levels[si + 1]:
                        substrate_vols[si].append(vtag)
                        break
                else:
                    air_vols.append(vtag)
            for si, sv in substrate_vols.items():
                if not sv:
                    raise RuntimeError(f"no volume found for substrate {si}")
                gmsh.model.addPhysicalGroup(
                    3, sv, self._TAG_SUBSTRATE_BASE + si)
            gmsh.model.addPhysicalGroup(3, air_vols, self._TAG_AIR)

            self._set_mesh_size(gmsh, lc, density, a, b, z_probe, min_t,
                                [t for sfs in patch_surfs.values()
                                 for t in sfs])

            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(msh_path)
        finally:
            gmsh.finalize()

    def _set_mesh_size(self, gmsh, lc, density, a, b, z_probe, min_t,
                       patch_tags):
        """Background size field: probe gap, patches, lc elsewhere.

        - Probe: size min(a, b - a) / 3 (3 elements across the narrower of
          the probe half-side and the gap) within the hole, growing with
          slope 0.5 away from the probe and hole edges.
        - Patches: size s = _PATCH_SIZE_X_DENSITY * min_t / density (h/3
          at density 12) inside the patch footprint from the ground plane
          to the highest patch, and within min_t of the patch outlines,
          growing with slope 0.5. s sets the accuracy more than lc: on
          the 2.45 GHz FR-4 patch (tools/patch_convergence.py, 10 MHz steps)
          f_res is 2.313, 2.360, 2.370 GHz at density 8, 12, 16 with
          s = 0.75 h, h/2, 0.375 h, and 2.377 GHz at density 12 with
          s = h/3, 2.388 GHz with s = h/4 (456k edges, 2.6x the time of
          h/3); openEMS gives 2.389 GHz.
        """
        gp = self.ground_plane_size
        cx = gp / 2 + self.probe_offset[0]
        cy = gp / 2 + self.probe_offset[1]
        tol = 1e-3 * min(a, min_t)

        probe_size = min(a, b - a) / 3.0
        probe_curves = [
            t for _, t in gmsh.model.getEntitiesInBoundingBox(
                cx - b - tol, cy - b - tol, -tol,
                cx + b + tol, cy + b + tol, z_probe + tol, 1)]
        patch_size = min(lc, self._PATCH_SIZE_X_DENSITY * min_t / density)
        patch_curves = [
            t for _, t in gmsh.model.getBoundary(
                [(2, t) for t in patch_tags], combined=True, oriented=False)]

        gmsh.option.setNumber("Mesh.MeshSizeMin", min(probe_size, patch_size))
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)

        field = gmsh.model.mesh.field
        field.add("Distance", 1)
        field.setNumbers(1, "CurvesList", probe_curves)
        field.setNumber(1, "Sampling", 20)
        field.add("Threshold", 2)
        field.setNumber(2, "InField", 1)
        field.setNumber(2, "SizeMin", probe_size)
        field.setNumber(2, "SizeMax", lc)
        field.setNumber(2, "DistMin", b - a)
        field.setNumber(2, "DistMax",
                        (b - a) + (lc - probe_size) / self._PROBE_SIZE_SLOPE)

        field.add("Distance", 3)
        field.setNumbers(3, "CurvesList", patch_curves)
        field.setNumber(3, "Sampling", 200)
        field.add("Threshold", 4)
        field.setNumber(4, "InField", 3)
        field.setNumber(4, "SizeMin", patch_size)
        field.setNumber(4, "SizeMax", lc)
        field.setNumber(4, "DistMin", min_t)
        field.setNumber(4, "DistMax",
                        min_t + (lc - patch_size) / self._PATCH_SIZE_SLOPE)

        x0 = min(gp / 2 - p['width'] / 2 for p in self.patches)
        y0 = min(gp / 2 - p['length'] / 2 for p in self.patches)
        field.add("Box", 6)
        field.setNumber(6, "VIn", patch_size)
        field.setNumber(6, "VOut", lc)
        field.setNumber(6, "XMin", x0)
        field.setNumber(6, "XMax", gp - x0)
        field.setNumber(6, "YMin", y0)
        field.setNumber(6, "YMax", gp - y0)
        field.setNumber(6, "ZMin", 0.0)
        field.setNumber(6, "ZMax", max(p['height'] for p in self.patches))

        field.add("Min", 5)
        field.setNumbers(5, "FieldsList", [2, 4, 6])
        field.setAsBackgroundMesh(5)

    def _setup_simulation(self, freq: float):
        """Set up boundary conditions and ports."""
        if self._mesh is None:
            raise RuntimeError("Mesh not generated. Call generate_mesh() first.")

        # Build PEC BC for all conducting surfaces
        bc = em.BC()
        pec_tags = [self._TAG_GROUND, self._TAG_PROBE]
        # Add patch PEC tags
        for pi in range(len(self.patches)):
            pec_tags.append(self._TAG_PATCH_BASE + pi)

        for tag in pec_tags:
            if em.has_surface_tag(self._mesh, tag):
                bc_tag = em.build_edge_pec(self._mesh, tag)
                bc.merge(bc_tag)

        self._bc = bc

        # Lumped port on the strip from the probe face to the hole edge at
        # z = 0; E along x, across the gap
        config = em.LumpedPortConfig()
        config.surface_tag = self._TAG_PORT
        config.z0 = self.z0
        config.e_direction = np.array([1.0, 0.0, 0.0])

        port = em.build_lumped_port(self._mesh, config)
        self._ports = [port]

    def simulate(self, freq: float) -> Tuple[complex, object]:
        """
        Run full-wave FEM simulation at given frequency.

        Uses calculate_sparams_with_fields() with the gap-port formulation
        of the C++ core (load j k0 eta0 / Z0 on the port strip); the lumped
        port needs no weight normalization. S11 and the field come from one
        factorization.

        Args:
            freq: Frequency in Hz

        Returns:
            Tuple of (S11, solution_vector). If the solve failed, S11 is NaN,
            the solution vector is empty and nothing is stored.
        """
        self._setup_simulation(freq)
        self._check_air_box(freq)
        omega = 2 * np.pi * freq

        params = em.MaxwellParams()
        params.omega = omega
        params.use_abc = True
        params.abc_surface_tags = {self._TAG_ABC}
        for si, sub in enumerate(self.substrates):
            tag = self._TAG_SUBSTRATE_BASE + si
            eps_r = sub['eps_r']
            tan_d = sub.get('tan_d', 0.0)
            params.set_eps_r_region(tag, complex(eps_r, -eps_r * tan_d))

        opts = em.SolveOptions()
        opts.use_direct = True

        # S-parameters and the port-0 field from one factorization
        S, fields = em.calculate_sparams_with_fields(
            self._mesh, params, self._bc, self._ports, opts
        )
        S11 = S[0, 0]
        solution = fields[0]  # empty when the solve failed (S11 is NaN)

        # Passivity diagnostic (warn but don't clamp)
        if abs(S11) > 1.01:
            import sys
            print(
                f"WARNING: |S11| = {abs(S11):.4f} > 1 at {freq/1e9:.3f} GHz "
                f"(passivity violation — consider refining mesh)",
                file=sys.stderr,
            )

        if len(solution) > 0:
            self._store_solution(freq, solution)
        return S11, solution

    def input_impedance(self, freq: float) -> complex:
        """
        Compute input impedance at given frequency.

        Args:
            freq: Frequency in Hz

        Returns:
            Complex input impedance in ohms
        """
        S11, _ = self.simulate(freq)
        Z0 = self.z0
        # Z_in = Z0 * (1 + S11) / (1 - S11)
        if abs(1 - S11) < 1e-30:
            return complex(1e6, 0)  # Open circuit
        return Z0 * (1 + S11) / (1 - S11)

    def return_loss(self, freq: float, z0: float = 50.0) -> float:
        """
        Compute return loss |S11| in dB.

        Args:
            freq: Frequency in Hz
            z0: Reference impedance

        Returns:
            Return loss in dB (negative value)
        """
        S11, _ = self.simulate(freq)
        return 20 * np.log10(abs(S11) + 1e-30)

    def convergence_check(
        self,
        freq: float,
        density: float = 10.0,
        refinement_factor: float = 1.5,
        threshold: float = 0.02,
        verbose: bool = True,
    ) -> dict:
        """Check mesh convergence by solving at two densities.

        Solves the problem at `density` and `density * refinement_factor`,
        then compares S11. If the difference exceeds `threshold`, the mesh
        is not converged and the user should increase density.

        Args:
            freq: Frequency in Hz.
            density: Coarse mesh density (elements per wavelength).
            refinement_factor: Multiplier for fine mesh density.
            threshold: Maximum |ΔS11| for convergence.
            verbose: Print convergence report.

        Returns:
            dict with keys: 'converged', 'S11_coarse', 'S11_fine',
            'delta', 'density_coarse', 'density_fine'.
        """
        # Solve at coarse density
        self.generate_mesh(density=density)
        S11_coarse, _ = self.simulate(freq)

        # Solve at fine density
        fine_density = density * refinement_factor
        self.generate_mesh(density=fine_density)
        S11_fine, _ = self.simulate(freq)

        delta = abs(S11_fine - S11_coarse)
        converged = delta < threshold

        result = {
            'converged': converged,
            'S11_coarse': S11_coarse,
            'S11_fine': S11_fine,
            'delta': delta,
            'density_coarse': density,
            'density_fine': fine_density,
        }

        if verbose:
            status = "CONVERGED" if converged else "NOT CONVERGED"
            print(f"Mesh convergence check at {freq/1e9:.3f} GHz: {status}")
            print(f"  Coarse (density={density:.0f}): "
                  f"|S11| = {abs(S11_coarse):.4f}")
            print(f"  Fine   (density={fine_density:.0f}): "
                  f"|S11| = {abs(S11_fine):.4f}")
            print(f"  ΔS11 = {delta:.4f} "
                  f"({'<' if converged else '>'} {threshold} threshold)")

        return result

    def frequency_sweep(
        self,
        f_start: float,
        f_stop: float,
        n_points: int = 21,
        verbose: bool = False,
        keep_fields: Union[None, str, float,
                           Sequence[Union[float, str]]] = None,
    ) -> FrequencySweepResult:
        """
        Sweep frequency and compute S11 and input impedance.

        Args:
            f_start: Start frequency (Hz)
            f_stop: Stop frequency (Hz)
            n_points: Number of frequency points
            verbose: Print progress
            keep_fields: Field solutions to keep for radiation_pattern(),
                directivity(), near_field() and export_fields_vtk(), which
                then need no solve at those frequencies. None (default)
                keeps none. "min_s11" keeps the field at the smallest
                |S11| (the first on a tie, which is where np.argmin of the
                returned |S11| lands); the sweep holds only the best field
                so far. A frequency or a sequence of frequencies (Hz)
                keeps those; each must be within 1e-9 relative of a sweep
                frequency. The sequence may also contain "min_s11". A kept
                field is stored under the sweep frequency from the returned
                array, so look it up with that value. It replaces any field
                already stored at that frequency.

        Returns:
            FrequencySweepResult: unpacks as (frequencies, Z_in_array,
            S11_array). With "min_s11" in keep_fields, its min_s11_index
            and min_s11_frequency name the frequency whose field was kept
            (None otherwise). Where a solve failed, S11 and Z_in are NaN
            and a RuntimeWarning names those frequencies.
        """
        freqs = np.linspace(f_start, f_stop, n_points)
        keep_indices, keep_min = self._keep_request(freqs, keep_fields)

        # Set up BC and ports (frequency-independent for the lumped port)
        f_mid = (f_start + f_stop) / 2.0
        self._setup_simulation(f_mid)
        self._check_air_box(float(freqs.min()))

        params = em.MaxwellParams()
        params.omega = 2 * np.pi * f_mid  # overwritten per frequency by the sweep
        params.use_abc = True
        params.abc_surface_tags = {self._TAG_ABC}
        for si, sub in enumerate(self.substrates):
            tag = self._TAG_SUBSTRATE_BASE + si
            eps_r = sub['eps_r']
            tan_d = sub.get('tan_d', 0.0)
            params.set_eps_r_region(tag, complex(eps_r, -eps_r * tan_d))

        opts = em.SolveOptions()
        opts.use_direct = True

        if verbose:
            print(f"  Running {n_points}-point frequency sweep "
                  f"({f_start/1e9:.3f}–{f_stop/1e9:.3f} GHz)")

        # Use efficient frequency sweep (K/M separation + symbolic reuse)
        sweep_result = em.frequency_sweep(
            self._mesh, params, self._bc, self._ports,
            freqs.tolist(), opts, keep_indices, keep_min
        )
        for kept in sweep_result.fields:
            x = kept.x[0]
            if len(x) > 0:  # empty when the solve at that frequency failed
                self._store_solution(float(freqs[kept.index]), x)

        Z_list = []
        S11_list = []
        failed = []
        Z0 = self.z0
        for i, S_mat in enumerate(sweep_result.S_matrices):
            S11 = complex(S_mat[0, 0])
            S11_list.append(S11)
            if np.isnan(S11.real) or np.isnan(S11.imag):
                # Failed solve: no S11, so no Z either.
                failed.append(float(freqs[i]))
                Z_in = complex(np.nan, np.nan)
            elif abs(1 - S11) > 1e-30:
                Z_in = Z0 * (1 + S11) / (1 - S11)
            else:
                Z_in = complex(1e6, 0)
            Z_list.append(Z_in)
        if failed:
            warnings.warn(
                f"StackedPatchDesign.frequency_sweep: the solve failed at "
                f"{len(failed)} of {len(freqs)} frequencies (S11 and Z_in "
                f"are NaN there): "
                + ", ".join(f"{f / 1e9:.6g} GHz" for f in failed),
                RuntimeWarning,
                stacklevel=2,
            )

        return FrequencySweepResult(freqs, np.array(Z_list),
                                    np.array(S11_list),
                                    sweep_result.min_reflection_index)

    @staticmethod
    def _keep_request(freqs: np.ndarray, keep_fields) -> Tuple[List[int], bool]:
        """Translate frequency_sweep()'s keep_fields into sweep indices and
        the min-|S11| flag. Raises ValueError before any solve."""
        if keep_fields is None:
            return [], False
        items = [keep_fields] if np.ndim(keep_fields) == 0 else keep_fields
        indices, keep_min = [], False
        for item in items:
            if isinstance(item, str):
                if item != "min_s11":
                    raise ValueError(
                        f"keep_fields: unknown option {item!r}; use "
                        f"'min_s11' or frequencies in Hz")
                keep_min = True
                continue
            f = float(item)
            i = int(np.argmin(np.abs(freqs - f)))
            if abs(freqs[i] - f) > 1e-9 * abs(f):
                raise ValueError(
                    f"keep_fields: {f:.6g} Hz is not a sweep frequency")
            indices.append(i)
        return indices, keep_min

    # Far-field evaluation works on blocks of directions so the
    # directions x faces phase matrix stays small. Per matrix entry a block
    # holds the real phase argument (8 bytes) and the complex exponential
    # (16 bytes); 32 allows for the temporaries of the matrix product.
    _FARFIELD_BYTES_PER_ENTRY = 32
    _FARFIELD_CHUNK_BYTES = 256 * 2**20

    def radiation_pattern(
        self,
        freq: float,
        n_theta: int = 91,
        n_phi: int = 73,
    ) -> "em.FFPattern3D":
        """
        Compute 3D radiation pattern via Huygens surface + Stratton-Chu.

        The pattern is cached per (freq, n_theta, n_phi) and reused until
        the field solution at freq is replaced or generate_mesh() is
        called; the cache holds the _PATTERN_CACHE_MAX most recently used
        patterns. Each call returns its own copy of the arrays.

        Args:
            freq: Frequency in Hz
            n_theta: Number of theta samples (0 to 90 deg)
            n_phi: Number of phi samples (0 to 360 deg)

        Returns:
            Pattern object with theta_grid, phi_grid, E_theta, E_phi
        """
        key = (freq, n_theta, n_phi)
        solution = self._solution_at(freq)
        cached = self._pattern_cache.pop(key, None)
        if cached is not None and cached[0] is solution:
            self._pattern_cache[key] = cached  # now most recently used
            return self._copy_pattern(cached[1])

        pattern = self._evaluate_pattern(freq, n_theta, n_phi)
        self._pattern_cache[key] = (solution, pattern)
        while len(self._pattern_cache) > self._PATTERN_CACHE_MAX:
            del self._pattern_cache[next(iter(self._pattern_cache))]
        return self._copy_pattern(pattern)

    @staticmethod
    def _copy_pattern(pattern):
        return _AperturePattern(
            pattern.theta_grid.copy(), pattern.phi_grid.copy(),
            pattern.E_theta.copy(), pattern.E_phi.copy(),
        )

    def _evaluate_pattern(self, freq: float, n_theta: int, n_phi: int):
        """Stratton-Chu sum over the Huygens faces for an n_theta x n_phi
        grid, vectorized over directions and faces in blocks."""
        solution = self._solutions[freq]
        omega = 2 * np.pi * freq
        c0 = 299792458.0
        k0 = omega / c0

        # Extract Huygens surface. Each attribute is read once: pybind
        # copies the whole face list on every access.
        huygens = em.extract_huygens_surface(
            self._mesh, solution, self._TAG_HUYGENS, omega
        )
        r = np.asarray(huygens.r, dtype=float).reshape(-1, 3)
        n = np.asarray(huygens.n, dtype=float).reshape(-1, 3)
        E_tan = np.asarray(huygens.E_tan, dtype=complex).reshape(-1, 3)
        area = np.asarray(huygens.area, dtype=float).reshape(-1)

        # Ensure consistent outward (+z) normals on the horizontal Huygens
        # surface. The C++ extractor picks one parent tet per face, and for
        # faces shared between tets above and below the surface the choice
        # is non-deterministic, occasionally flipping the normal.
        # Use aperture/image theory: for an antenna above a PEC ground
        # plane, the far field equals 2x the magnetic-current (M = -n x E)
        # contribution from the aperture. H_tan is not used because image
        # theory accounts for it via the factor of 2.
        sign = np.where(n[:, 2] < 0.0, -1.0, 1.0)[:, None]
        M = -np.cross(sign * n, sign * 2.0 * E_tan)
        weighted_M = M * area[:, None]                       # (faces, 3)

        # Upper hemisphere only (ground plane blocks lower hemisphere)
        theta_rad = np.linspace(0, np.pi / 2, n_theta)
        phi_rad = np.linspace(0, 2 * np.pi, n_phi)
        theta_grid, phi_grid = np.meshgrid(theta_rad, phi_rad, indexing="ij")
        th = theta_grid.ravel()
        ph = phi_grid.ravel()
        sin_th, cos_th = np.sin(th), np.cos(th)
        sin_ph, cos_ph = np.sin(ph), np.cos(ph)
        zeros = np.zeros_like(th)
        rhat = np.stack([sin_th * cos_ph, sin_th * sin_ph, cos_th], axis=1)
        th_hat = np.stack([cos_th * cos_ph, cos_th * sin_ph, -sin_th], axis=1)
        ph_hat = np.stack([-sin_ph, cos_ph, zeros], axis=1)

        # Per direction, the loop form summed
        #   j k0 (rhat x M_s) x rhat * exp(-j k0 rhat.r_s) * area_s * cos_th
        # over faces s. (rhat x M) x rhat is the part of M perpendicular to
        # rhat, and th_hat, ph_hat are perpendicular to rhat, so
        #   E_theta = j k0 cos_th * th_hat . F,  E_phi = j k0 cos_th * ph_hat . F
        # with F = sum_s exp(-j k0 rhat.r_s) area_s M_s. The aperture
        # obliquity factor cos_th approximates the horizon cancellation a
        # closed Huygens surface would give; an open single-plane surface
        # has no side faces.
        n_dir = th.size
        n_faces = r.shape[0]
        per_dir = self._FARFIELD_BYTES_PER_ENTRY * max(n_faces, 1)
        chunk = max(1, int(self._FARFIELD_CHUNK_BYTES // per_dir))

        F = np.zeros((n_dir, 3), dtype=complex)
        if n_faces:
            for start in range(0, n_dir, chunk):
                stop = min(start + chunk, n_dir)
                arg = rhat[start:stop] @ r.T                  # (dirs, faces)
                arg *= -k0
                phase = arg * 1j                              # complex copy
                del arg
                np.exp(phase, out=phase)
                F[start:stop] = phase @ weighted_M

        scale = 1j * k0 * cos_th
        E_theta_out = (scale * np.einsum("ij,ij->i", F, th_hat)).reshape(
            n_theta, n_phi)
        E_phi_out = (scale * np.einsum("ij,ij->i", F, ph_hat)).reshape(
            n_theta, n_phi)

        # Lightweight pattern object compatible with the C++ API
        return _AperturePattern(theta_grid, phi_grid, E_theta_out, E_phi_out)

    def directivity(
        self,
        freq: float,
        n_theta: int = 91,
        n_phi: int = 73,
    ) -> float:
        """
        Compute directivity at given frequency.

        The value depends on the pattern grid only through the numerical
        integration, and the default grid is fixed (91 x 73), so it does
        not depend on which patterns were computed earlier. The pattern
        comes from the same cache as radiation_pattern(): a call with the
        grid of an earlier radiation_pattern() call does not evaluate the
        far field again.

        Args:
            freq: Frequency in Hz
            n_theta, n_phi: Pattern grid used for the integration

        Returns:
            Directivity (linear scale)
        """
        pattern = self.radiation_pattern(freq, n_theta=n_theta, n_phi=n_phi)
        if isinstance(pattern, _AperturePattern):
            return self._compute_directivity_py(pattern)
        return em.compute_directivity(pattern)

    @staticmethod
    def _compute_directivity_py(pattern):
        """Compute directivity from a Python pattern (upper hemisphere)."""
        Z0 = 376.73
        pwr = pattern.power_pattern()
        U_max = pwr.max() / (2.0 * Z0)
        if U_max < 1e-30:
            return 0.0
        Nth, Nph = pwr.shape
        if Nth < 2 or Nph < 2:
            return 1.0
        dtheta = (pattern.theta_grid[-1, 0] - pattern.theta_grid[0, 0]) / (Nth - 1)
        dphi = (pattern.phi_grid[0, -1] - pattern.phi_grid[0, 0]) / (Nph - 1)
        P_rad = 0.0
        for i in range(Nth):
            theta = pattern.theta_grid[i, 0]
            sin_th = np.sin(theta)
            wt_th = 0.5 if (i == 0 or i == Nth - 1) else 1.0
            for j in range(Nph):
                wt_ph = 0.5 if (j == 0 or j == Nph - 1) else 1.0
                U = pwr[i, j] / (2.0 * Z0)
                P_rad += U * sin_th * wt_th * wt_ph
        P_rad *= dtheta * dphi
        # For upper hemisphere pattern, assume lower hemisphere has zero
        # contribution (ground plane). Multiply by 2 for 4*pi normalization.
        if P_rad < 1e-30:
            return 0.0
        return 4.0 * np.pi * U_max / (2.0 * P_rad)

    def near_field(
        self, freq: float, points: np.ndarray
    ) -> np.ndarray:
        """
        Evaluate E-field at specified points.

        Args:
            freq: Frequency in Hz
            points: Nx3 array of evaluation points

        Returns:
            Nx3 complex array of E-field values
        """
        solution = self._solution_at(freq)
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        E_out = np.zeros((len(points), 3), dtype=complex)

        # Mesh arrays, read once (the list attributes copy on every access)
        mesh = self._mesh
        xyz = mesh.node_coords()
        tet_verts = xyz[mesh.tet_nodes()]                    # (tets, 4, 3)
        tet_edges = mesh.tet_edges()
        tet_orient = mesh.tet_edge_orient()

        # Barycentric coordinates of each point in every tet, in numpy, pick
        # the candidates; the first candidate in mesh order that passes the
        # C++ test (all coordinates >= -1e-10) holds the point, as in the
        # former tet-by-tet search. The numpy prefilter is looser (-1e-6) so
        # round-off cannot drop the tet the C++ test would accept.
        T = np.stack([tet_verts[:, k] - tet_verts[:, 3] for k in range(3)],
                     axis=2)                                  # (tets, 3, 3)
        T_inv = np.linalg.pinv(T)   # pinv: a flat tet does not raise
        for pi, point in enumerate(points):
            lam012 = np.einsum('tij,tj->ti', T_inv, point - tet_verts[:, 3])
            lam_min = np.minimum(lam012.min(axis=1),
                                 1.0 - lam012.sum(axis=1))
            for t in np.flatnonzero(lam_min >= -1e-6):
                verts = [tet_verts[t, k] for k in range(4)]
                lam = em.compute_barycentric(verts, point)
                if all(v >= -1e-10 for v in lam):
                    edge_dofs = solution[tet_edges[t].tolist()]
                    E_out[pi] = em.evaluate_edge_field(
                        verts, tet_orient[t].tolist(), edge_dofs, point)
                    break

        return E_out

    def export_fields_vtk(self, freq: float, path: str) -> None:
        """
        Export E-field solution to VTK file.

        Args:
            freq: Frequency in Hz
            path: Output file path (.vtu)
        """
        em.export_fields_vtk(self._mesh, self._solution_at(freq), path)

    def export_touchstone(
        self, path: str, freqs: np.ndarray, S11_array: np.ndarray
    ) -> None:
        """
        Export S11 to Touchstone .s1p format.

        Args:
            path: Output path
            freqs: Frequencies in Hz
            S11_array: Complex S11 values
        """
        S_matrices = [np.array([[s]]) for s in S11_array]
        if not path.endswith(".s1p"):
            path = path + ".s1p"
        em.write_touchstone_nport(path, list(freqs), S_matrices)

    def __repr__(self) -> str:
        n_patches = len(self.patches)
        n_subs = len(self.substrates)
        return (
            f"StackedPatchDesign({n_patches} patch{'es' if n_patches > 1 else ''}, "
            f"{n_subs} substrate{'s' if n_subs > 1 else ''}, "
            f"probe r={self.probe_radius*1000:.2f}mm, "
            f"z0={self.z0:.0f}ohm)"
        )
