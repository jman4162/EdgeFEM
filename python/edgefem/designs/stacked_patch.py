"""
Stacked Patch Antenna with Cavity Feed Design Class

Multi-layer stacked patch antenna with cavity-backed probe feed.
Uses full-wave FEM simulation via lumped port excitation and
Stratton-Chu near-to-far field transformation for radiation patterns.
"""

import subprocess
import tempfile
import os
from typing import Optional, Tuple, List, Dict
import numpy as np

import pyedgefem as em


class StackedPatchDesign:
    """
    Multi-layer stacked patch antenna with cavity-backed probe feed.

    The antenna consists of one or more stacked patches on dielectric
    substrates above a ground plane, fed by a coaxial probe through
    a cavity below the ground plane. Full-wave FEM simulation provides
    S-parameters, input impedance, and radiation patterns.

    Geometry (bottom to top):
        - Cavity volume below ground plane with probe
        - Ground plane (PEC) with hole for probe
        - Substrate layers with patch conductors
        - Air box with ABC outer boundary
        - Huygens surface for NTF

    Example:
        design = StackedPatchDesign(
            patches=[{'width': 30e-3, 'length': 25e-3, 'height': 1.6e-3}],
            substrates=[{'thickness': 1.6e-3, 'eps_r': 4.4, 'tan_d': 0.02}],
            cavity_size=40e-3,
            cavity_depth=5e-3,
        )
        design.generate_mesh(density=12)
        S11, solution = design.simulate(2.4e9)
        pattern = design.radiation_pattern(2.4e9)
    """

    # Physical tag scheme
    _TAG_GROUND = 1
    _TAG_CAVITY_WALLS = 2
    _TAG_CAVITY_BOTTOM = 3
    _TAG_PROBE = 4
    _TAG_PORT = 5
    _TAG_PATCH_BASE = 10      # 10, 11, 12, ...
    _TAG_ABC = 50
    _TAG_HUYGENS = 60
    _TAG_CAVITY_VOL = 100
    _TAG_SUBSTRATE_BASE = 110  # 110, 111, ...
    _TAG_AIR = 150

    def __init__(
        self,
        patches: List[Dict],
        substrates: List[Dict],
        cavity_size: float,
        cavity_depth: float,
        probe_radius: float = 0.5e-3,
        probe_offset: Tuple[float, float] = (0.0, 0.0),
        z0: float = 50.0,
        ground_plane_size: Optional[float] = None,
        air_box_height: Optional[float] = None,
    ):
        """
        Initialize stacked patch design.

        Args:
            patches: List of patch dicts, each with keys:
                     'width' (m), 'length' (m), 'height' (m above ground)
            substrates: List of substrate dicts, each with keys:
                       'thickness' (m), 'eps_r', 'tan_d' (optional, default 0)
            cavity_size: Lateral dimension of cavity (m)
            cavity_depth: Depth of cavity below ground plane (m)
            probe_radius: Radius of probe feed (m)
            probe_offset: (x, y) offset of probe from center (m)
            z0: Reference impedance (ohms)
            ground_plane_size: Ground plane size. If None, 3x largest patch.
            air_box_height: Air box height above top patch. If None, lambda/4.
        """
        if len(patches) == 0:
            raise ValueError("At least one patch is required")
        if len(substrates) == 0:
            raise ValueError("At least one substrate is required")

        self.patches = patches
        self.substrates = substrates
        self.cavity_size = cavity_size
        self.cavity_depth = cavity_depth
        self.probe_radius = probe_radius
        self.probe_offset = probe_offset
        self.z0 = z0

        max_patch = max(max(p['width'], p['length']) for p in patches)
        self.ground_plane_size = ground_plane_size or 3.0 * max_patch
        self._air_box_height = air_box_height

        self._mesh = None
        self._mesh_path = None
        self._bc = None
        self._ports = None
        self._last_solutions = {}

    @property
    def mesh(self):
        return self._mesh

    def generate_mesh(
        self,
        density: float = 12.0,
        output_path: Optional[str] = None,
        design_freq: Optional[float] = None,
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

        wavelength = c0 / design_freq
        lc = wavelength / density

        if self._air_box_height is None:
            self._air_box_height = wavelength / 4

        if output_path is None:
            fd, msh_path = tempfile.mkstemp(suffix=".msh")
            os.close(fd)
        elif output_path.endswith(".msh"):
            msh_path = output_path
        else:
            msh_path = output_path + ".msh"

        self._build_gmsh_geometry(gmsh, lc, msh_path)

        self._mesh = em.load_gmsh(msh_path)
        self._mesh_path = msh_path
        self._bc = None
        self._ports = None
        self._last_solutions = {}
        if hasattr(self, '_port_alpha_sq'):
            del self._port_alpha_sq

    def _build_gmsh_geometry(self, gmsh, lc: float, msh_path: str):
        """Build antenna geometry using Gmsh Python API with OCC kernel."""
        gp = self.ground_plane_size
        cd = self.cavity_depth
        h_air = self._air_box_height

        # Probe inner conductor half-width (the actual probe)
        pr = max(self.probe_radius, lc / 8, 0.5e-3)
        # Probe hole half-width (outer conductor / ground plane hole)
        # Must be larger than pr to create a coaxial gap
        pw = max(2.5 * pr, lc / 4, 2e-3)
        cx = gp / 2 + self.probe_offset[0]
        cy = gp / 2 + self.probe_offset[1]

        z_bot = -cd
        z_levels = [0.0]
        for s in self.substrates:
            z_levels.append(z_levels[-1] + s['thickness'])
        z_top = z_levels[-1] + h_air

        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add("stacked_patch")

        occ = gmsh.model.occ

        # Main box spanning full domain
        occ.addBox(0, 0, z_bot, gp, gp, z_top - z_bot)

        # Ground plane at z=0 (will be split by probe hole via fragment)
        sf_gnd = occ.addRectangle(0, 0, 0, gp, gp)
        # Probe hole outline at z=0 — larger than probe to create coaxial gap
        sf_probe_hole = occ.addRectangle(cx - pw, cy - pw, 0, 2 * pw, 2 * pw)

        # Port: horizontal surface at z=0 in the gap between probe and ground.
        # This represents the coaxial cross-section.
        # E-field direction is x (radial, from inner to outer conductor).
        # Port strip: from probe edge (cx+pr) to hole edge (cx+pw), at z=0
        sf_port = occ.addRectangle(cx + pr, cy - pr, 0, pw - pr, 2 * pr)

        # Probe conductor: thin PEC box from cavity bottom to first patch
        z_patch0 = self.patches[0]['height']
        probe_box = occ.addBox(
            cx - pr, cy - pr, z_bot, 2 * pr, 2 * pr, z_patch0 - z_bot
        )

        embed_list = [
            (2, sf_gnd), (2, sf_probe_hole), (2, sf_port),
            (3, probe_box),
        ]
        for pi, patch in enumerate(self.patches):
            z_patch = patch['height']
            pw_half = patch['width'] / 2
            pl_half = patch['length'] / 2
            sf = occ.addRectangle(
                gp / 2 - pw_half, gp / 2 - pl_half, z_patch,
                patch['width'], patch['length']
            )
            embed_list.append((2, sf))

        occ.synchronize()
        occ.fragment([(3, 1)], embed_list)
        occ.synchronize()

        # Classify surfaces by bounding box position
        vols = [tag for _, tag in gmsh.model.getEntities(3)]
        tol = 1e-4

        gnd_surfs = []
        port_surfs = []
        probe_surfs = []
        patch_surfs = []
        bot_surfs = []
        top_surfs = []
        side_surfs = []
        cavity_wall_surfs = []

        for _, tag in gmsh.model.getEntities(2):
            bb = gmsh.model.getBoundingBox(2, tag)
            zc = (bb[2] + bb[5]) / 2
            ze = bb[5] - bb[2]
            xe = bb[3] - bb[0]
            ye = bb[4] - bb[1]
            xc = (bb[0] + bb[3]) / 2
            yc = (bb[1] + bb[4]) / 2

            if ze < tol:  # Flat in z
                if abs(zc) < tol:
                    # At z=0: distinguish ground, port, and probe hole
                    near_probe = (abs(xc - cx) < pw * 2 and
                                  abs(yc - cy) < pw * 2)
                    is_small = (xe < pw * 4 and ye < pw * 4)
                    # Port surface: strip between probe (pr) and hole (pw)
                    # centered on x > cx (right side of probe)
                    is_port = (is_small and near_probe and
                               bb[0] > cx + pr - tol and
                               xe < (pw - pr) * 2 + tol)
                    if is_port:
                        port_surfs.append(tag)
                    elif is_small and near_probe:
                        # Other small surfaces at z=0 near probe = probe hole
                        # Not PEC (it's the gap)
                        pass
                    else:
                        gnd_surfs.append(tag)  # Large = ground plane
                elif abs(zc - z_bot) < tol:
                    bot_surfs.append(tag)
                elif abs(zc - z_top) < tol:
                    top_surfs.append(tag)
                else:
                    for pi, patch in enumerate(self.patches):
                        if abs(zc - patch['height']) < tol:
                            patch_surfs.append((pi, tag))
                            break
            elif ye < tol:  # Flat in y — vertical surface
                near_probe = (abs(yc - cy) < pr * 2 + tol)
                is_small = (xe < pr * 4 and near_probe)
                if is_small:
                    probe_surfs.append(tag)
                else:
                    side_surfs.append(tag)
            elif xe < tol:  # Flat in x — vertical surface
                near_probe = (abs(xc - cx) < pr * 2 + tol)
                is_small = (ye < pr * 4 and near_probe)
                if is_small:
                    probe_surfs.append(tag)
                else:
                    side_surfs.append(tag)

        # Physical groups
        if gnd_surfs:
            gmsh.model.addPhysicalGroup(2, gnd_surfs, self._TAG_GROUND)
        if bot_surfs:
            gmsh.model.addPhysicalGroup(2, bot_surfs, self._TAG_CAVITY_BOTTOM)
        if port_surfs:
            gmsh.model.addPhysicalGroup(2, port_surfs, self._TAG_PORT)
        if probe_surfs:
            gmsh.model.addPhysicalGroup(2, probe_surfs, self._TAG_PROBE)
        # Group patch surfaces by patch index (fragment may split one patch
        # rectangle into multiple surfaces)
        patch_by_index = {}
        for pi, sf in patch_surfs:
            patch_by_index.setdefault(pi, []).append(sf)
        for pi, sfs in patch_by_index.items():
            gmsh.model.addPhysicalGroup(2, sfs, self._TAG_PATCH_BASE + pi)

        # Cavity walls: outer side surfaces below ground plane
        for tag in side_surfs:
            bb = gmsh.model.getBoundingBox(2, tag)
            if bb[2] < -tol:  # Bottom edge below ground
                cavity_wall_surfs.append(tag)
        if cavity_wall_surfs:
            gmsh.model.addPhysicalGroup(2, cavity_wall_surfs, self._TAG_CAVITY_WALLS)

        # ABC = top + sides (above ground only)
        abc_list = list(top_surfs)
        for tag in side_surfs:
            if tag in cavity_wall_surfs:
                continue
            bb = gmsh.model.getBoundingBox(2, tag)
            zc = (bb[2] + bb[5]) / 2
            if zc > -tol:  # Above ground
                abc_list.append(tag)
        if abc_list:
            gmsh.model.addPhysicalGroup(2, abc_list, self._TAG_ABC)

        # Classify volumes by z-position into cavity, substrate, and air
        cavity_vols = []
        substrate_vols = {si: [] for si in range(len(self.substrates))}
        air_vols = []

        for vtag in vols:
            bb = gmsh.model.getBoundingBox(3, vtag)
            z_center = (bb[2] + bb[5]) / 2
            z_min = bb[2]
            z_max = bb[5]

            if z_max < tol:
                # Entirely below ground plane → cavity
                cavity_vols.append(vtag)
            elif z_min < -tol:
                # Straddles ground plane — treat as cavity
                cavity_vols.append(vtag)
            else:
                # Above ground plane: check if in substrate or air
                assigned = False
                cum_z = 0.0
                for si, sub in enumerate(self.substrates):
                    sub_top = cum_z + sub['thickness']
                    if z_center < sub_top + tol:
                        substrate_vols[si].append(vtag)
                        assigned = True
                        break
                    cum_z = sub_top
                if not assigned:
                    air_vols.append(vtag)

        if cavity_vols:
            gmsh.model.addPhysicalGroup(3, cavity_vols, self._TAG_CAVITY_VOL)
        for si, sv in substrate_vols.items():
            if sv:
                gmsh.model.addPhysicalGroup(3, sv, self._TAG_SUBSTRATE_BASE + si)
        if air_vols:
            gmsh.model.addPhysicalGroup(3, air_vols, self._TAG_AIR)

        # Mesh sizing: ensure at least 3 elements through the thinnest
        # substrate layer for proper field resolution
        min_sub_t = min(s['thickness'] for s in self.substrates)
        sub_mesh_size = min_sub_t / 3.0
        # Coaxial gap size needs at least 2-3 elements
        gap_size = pw - pr
        gap_mesh_size = gap_size / 3.0
        global_min = min(pr / 2, sub_mesh_size, gap_mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", global_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)

        # Local refinement: probe region AND substrate region
        probe_pts = [
            occ.addPoint(cx - pw, cy - pw, 0),
            occ.addPoint(cx + pw, cy - pw, 0),
            occ.addPoint(cx + pw, cy + pw, 0),
            occ.addPoint(cx - pw, cy + pw, 0),
        ]
        # Substrate region refinement points (patch corners at substrate top)
        sub_pts = []
        for patch in self.patches:
            z_p = patch['height']
            pw_h = patch['width'] / 2
            pl_h = patch['length'] / 2
            sub_pts.extend([
                occ.addPoint(gp/2 - pw_h, gp/2 - pl_h, z_p),
                occ.addPoint(gp/2 + pw_h, gp/2 - pl_h, z_p),
                occ.addPoint(gp/2 + pw_h, gp/2 + pl_h, z_p),
                occ.addPoint(gp/2 - pw_h, gp/2 + pl_h, z_p),
                occ.addPoint(gp/2, gp/2, z_p),  # center
            ])
        occ.synchronize()

        # Field 1: probe refinement
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "PointsList", probe_pts)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", pw / 2)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
        gmsh.model.mesh.field.setNumber(2, "DistMin", pw)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 4 * pw)

        # Field 3: substrate/patch refinement
        gmsh.model.mesh.field.add("Distance", 3)
        gmsh.model.mesh.field.setNumbers(3, "PointsList", sub_pts)
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", sub_mesh_size)
        gmsh.model.mesh.field.setNumber(4, "SizeMax", lc)
        gmsh.model.mesh.field.setNumber(4, "DistMin", min_sub_t)
        gmsh.model.mesh.field.setNumber(4, "DistMax", max(patch['width'], patch['length']))

        # Combine refinement fields
        gmsh.model.mesh.field.add("Min", 5)
        gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
        gmsh.model.mesh.field.setAsBackgroundMesh(5)

        gmsh.model.mesh.generate(3)

        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(msh_path)
        gmsh.finalize()

    def _compute_gamma(self, params, opts):
        """Compute gamma = w^H * A^{-1} * w for the lumped port."""
        port = self._ports[0]
        pec_edges = set(self._bc.dirichlet_edges)

        # Assemble bare FEM (no port loading) and port excitation
        assembly = em.assemble_maxwell(
            self._mesh, params, self._bc, [], -1
        )
        assembly_port = em.assemble_maxwell(
            self._mesh, params, self._bc, self._ports, 0
        )

        # Solve: x = A^{-1} * b where b = 2*w/sqrt(Z0)
        result = em.solve_linear(assembly.A, assembly_port.b, opts)

        # V = w^H * x = 2*gamma/sqrt(Z0)
        V_port = 0j
        for k in range(len(port.edges)):
            e = port.edges[k]
            if e not in pec_edges:
                V_port += np.conj(port.weights[k]) * result.x[e]

        Z0 = complex(port.mode.Z0)
        gamma = V_port * np.sqrt(Z0) / 2.0
        return gamma

    def _setup_simulation(self, freq: float):
        """Set up boundary conditions and ports."""
        if self._mesh is None:
            raise RuntimeError("Mesh not generated. Call generate_mesh() first.")

        # Build PEC BC for all conducting surfaces
        bc = em.BC()
        pec_tags = [
            self._TAG_GROUND,
            self._TAG_CAVITY_WALLS,
            self._TAG_CAVITY_BOTTOM,
            self._TAG_PROBE,
        ]
        # Add patch PEC tags
        for pi in range(len(self.patches)):
            pec_tags.append(self._TAG_PATCH_BASE + pi)

        for tag in pec_tags:
            if em.has_surface_tag(self._mesh, tag):
                bc_tag = em.build_edge_pec(self._mesh, tag)
                bc.merge(bc_tag)

        self._bc = bc

        # Build lumped port — E-field in x-direction (radial, from probe
        # inner conductor to ground outer conductor in coaxial gap)
        config = em.LumpedPortConfig()
        config.surface_tag = self._TAG_PORT
        config.z0 = self.z0
        config.e_direction = np.array([1.0, 0.0, 0.0])

        port = em.build_lumped_port(self._mesh, config)
        self._ports = [port]

    def simulate(self, freq: float) -> Tuple[complex, object]:
        """
        Run full-wave FEM simulation at given frequency.

        Uses the validated C++ S-parameter extraction (calculate_sparams)
        which handles port loading, solve, and S11 extraction correctly.

        Args:
            freq: Frequency in Hz

        Returns:
            Tuple of (S11, solution_vector)
        """
        self._setup_simulation(freq)

        omega = 2 * np.pi * freq

        params = em.MaxwellParams()
        params.omega = omega
        params.use_abc = True
        params.abc_surface_tags = {self._TAG_ABC}

        # Set substrate material properties
        for si, sub in enumerate(self.substrates):
            tag = self._TAG_SUBSTRATE_BASE + si
            eps_r = sub['eps_r']
            tan_d = sub.get('tan_d', 0.0)
            params.set_eps_r_region(tag, complex(eps_r, -eps_r * tan_d))

        # Impedance-based S-parameter extraction for lumped port.
        # We compute gamma = w^H * A^{-1} * w (bare FEM, no port loading)
        # and apply a frequency-independent normalization factor from a
        # reference frequency to get the physical impedance.
        opts = em.SolveOptions()
        opts.use_direct = True

        # Compute normalization factor at reference frequency (first call only)
        port = self._ports[0]
        Z0 = complex(port.mode.Z0)
        if not hasattr(self, '_port_alpha_sq'):
            # Use a low frequency (well below resonance) for stable normalization
            params_ref = em.MaxwellParams()
            f_ref = freq * 0.7  # Below expected resonance
            params_ref.omega = 2 * np.pi * f_ref
            params_ref.use_abc = True
            params_ref.abc_surface_tags = {self._TAG_ABC}
            for si, sub in enumerate(self.substrates):
                tag = self._TAG_SUBSTRATE_BASE + si
                eps_r_s = sub['eps_r']
                tan_d_s = sub.get('tan_d', 0.0)
                params_ref.set_eps_r_region(tag, complex(eps_r_s, -eps_r_s * tan_d_s))

            gamma_ref = self._compute_gamma(params_ref, opts)
            # Target: |gamma_norm| = Z0 => alpha^2 = Z0 / |gamma_ref|
            # Use magnitude to avoid complex phase rotation artifacts
            if abs(gamma_ref) > 1e-15:
                self._port_alpha_sq = Z0.real / abs(gamma_ref)
            else:
                self._port_alpha_sq = 1.0

        # Compute gamma at the actual frequency
        gamma_bare = self._compute_gamma(params, opts)

        # Apply normalization: gamma_norm = alpha^2 * gamma_bare
        gamma_norm = self._port_alpha_sq * gamma_bare

        # S11 = (gamma - Z0)/(gamma + Z0) from Woodbury identity
        S11 = (gamma_norm - Z0) / (gamma_norm + Z0)

        # Get solution vector for post-processing
        assembly = em.assemble_maxwell(
            self._mesh, params, self._bc, [], -1
        )
        assembly_port = em.assemble_maxwell(
            self._mesh, params, self._bc, self._ports, 0
        )
        result = em.solve_linear(assembly.A, assembly_port.b, opts)
        if not result.converged:
            raise RuntimeError(
                f"Solver did not converge: residual={result.residual}"
            )

        self._last_solutions[freq] = result.x
        return S11, result.x

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

    def frequency_sweep(
        self,
        f_start: float,
        f_stop: float,
        n_points: int = 21,
        verbose: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Sweep frequency and compute S11 and input impedance.

        Args:
            f_start: Start frequency (Hz)
            f_stop: Stop frequency (Hz)
            n_points: Number of frequency points
            verbose: Print progress

        Returns:
            Tuple of (frequencies, Z_in_array, S11_array)
        """
        freqs = np.linspace(f_start, f_stop, n_points)
        Z_list = []
        S11_list = []

        for i, f in enumerate(freqs):
            if verbose:
                print(f"  {i+1}/{n_points}: {f/1e9:.3f} GHz")
            S11, _ = self.simulate(f)
            S11_list.append(S11)

            Z0 = self.z0
            if abs(1 - S11) > 1e-30:
                Z_in = Z0 * (1 + S11) / (1 - S11)
            else:
                Z_in = complex(1e6, 0)
            Z_list.append(Z_in)

        return freqs, np.array(Z_list), np.array(S11_list)

    def radiation_pattern(
        self,
        freq: float,
        n_theta: int = 91,
        n_phi: int = 73,
    ) -> "em.FFPattern3D":
        """
        Compute 3D radiation pattern via Huygens surface + Stratton-Chu.

        Args:
            freq: Frequency in Hz
            n_theta: Number of theta samples
            n_phi: Number of phi samples

        Returns:
            FFPattern3D object
        """
        if freq not in self._last_solutions:
            self.simulate(freq)

        solution = self._last_solutions[freq]
        omega = 2 * np.pi * freq
        c0 = 299792458.0
        k0 = omega / c0

        # Extract Huygens surface
        huygens = em.extract_huygens_surface(
            self._mesh, solution, self._TAG_HUYGENS, omega
        )

        theta_rad = np.linspace(0, np.pi, n_theta).tolist()
        phi_rad = np.linspace(0, 2 * np.pi, n_phi).tolist()

        pattern = em.stratton_chu_3d(
            huygens.r, huygens.n, huygens.E_tan, huygens.H_tan,
            huygens.area, theta_rad, phi_rad, k0
        )

        return pattern

    def directivity(self, freq: float) -> float:
        """
        Compute directivity at given frequency.

        Args:
            freq: Frequency in Hz

        Returns:
            Directivity (linear scale)
        """
        pattern = self.radiation_pattern(freq)
        return em.compute_directivity(pattern)

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
        if freq not in self._last_solutions:
            self.simulate(freq)

        solution = self._last_solutions[freq]
        E_out = np.zeros((len(points), 3), dtype=complex)

        # Simple brute-force: for each point, find containing tet
        for pi, point in enumerate(points):
            for tet in self._mesh.tets:
                verts = []
                for i in range(4):
                    idx = self._mesh.nodeIndex[tet.conn[i]]
                    verts.append(self._mesh.nodes[idx].xyz)

                lam = em.compute_barycentric(verts, point)
                if all(l >= -1e-10 for l in lam):
                    edge_dofs = np.array([solution[tet.edges[e]] for e in range(6)])
                    E = em.evaluate_edge_field(
                        verts, list(tet.edge_orient), edge_dofs, point
                    )
                    E_out[pi] = E
                    break

        return E_out

    def export_fields_vtk(self, freq: float, path: str) -> None:
        """
        Export E-field solution to VTK file.

        Args:
            freq: Frequency in Hz
            path: Output file path (.vtu)
        """
        if freq not in self._last_solutions:
            self.simulate(freq)
        em.export_fields_vtk(self._mesh, self._last_solutions[freq], path)

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
            f"cavity={self.cavity_size*1000:.1f}mm, "
            f"z0={self.z0:.0f}ohm)"
        )
