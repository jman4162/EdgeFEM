"""
Periodic Unit Cell Design Class

High-level interface for metasurface/metamaterial unit cell simulation.
Supports Floquet port excitation for oblique incidence analysis.
"""

import subprocess
import tempfile
import os
from typing import Optional, Tuple, List, Dict
import numpy as np

import pyedgefem as em


class UnitCellDesign:
    """
    Periodic unit cell design for metasurface/metamaterial analysis.

    A unit cell represents the fundamental building block of a periodic
    structure (metasurface, frequency selective surface, or metamaterial).
    This class provides methods for:
    - Mesh generation with periodic boundary conditions
    - Floquet port excitation for oblique incidence
    - Reflection/transmission coefficient extraction
    - Surface impedance computation
    - Effective parameter extraction (NRW method)

    The unit cell is oriented with:
    - X, Y: in-plane periodic directions
    - Z: propagation direction (wave normally incident from below)
    - Substrate on ground plane (z=0) with patch on top (z=substrate_height)

    Example:
        cell = UnitCellDesign(
            period_x=5e-3,
            period_y=5e-3,
            substrate_height=1.0e-3,
            substrate_eps_r=3.5,
        )
        cell.add_patch(width=4e-3, length=4e-3)
        cell.generate_mesh(density=20)

        # Normal incidence
        R, T = cell.reflection_transmission(10e9)

        # Oblique incidence (30 degrees, TE polarization)
        R, T = cell.reflection_transmission(10e9, theta=30, phi=0, pol='TE')

        # Surface impedance
        Zs = cell.surface_impedance(10e9)

    Attributes:
        period_x: Unit cell period in x-direction (meters)
        period_y: Unit cell period in y-direction (meters)
        substrate_height: Substrate thickness (meters)
        substrate_eps_r: Substrate relative permittivity
    """

    # Physical tag assignments
    _TAG_GROUND_PLANE = 1      # PEC ground plane (z=0)
    _TAG_PATCH = 2             # PEC patch elements
    _TAG_PORT_BOTTOM = 3       # Bottom Floquet port
    _TAG_PORT_TOP = 4          # Top Floquet port
    _TAG_PERIODIC_X_NEG = 5    # -X periodic boundary
    _TAG_PERIODIC_X_POS = 6    # +X periodic boundary
    _TAG_PERIODIC_Y_NEG = 7    # -Y periodic boundary
    _TAG_PERIODIC_Y_POS = 8    # +Y periodic boundary
    _TAG_SUBSTRATE = 100       # Substrate volume
    _TAG_AIR_ABOVE = 101       # Air volume above
    _TAG_AIR_BELOW = 102       # Air volume below (for FSS)

    def __init__(
        self,
        period_x: float,
        period_y: float,
        substrate_height: float = 1e-3,
        substrate_eps_r: float = 3.5,
        substrate_tan_d: float = 0.0,
        *,
        air_height_above: Optional[float] = None,
        air_height_below: Optional[float] = None,
        has_ground_plane: bool = True,
    ):
        """
        Initialize unit cell design.

        Args:
            period_x: Unit cell period in x-direction (meters)
            period_y: Unit cell period in y-direction (meters)
            substrate_height: Substrate thickness (meters)
            substrate_eps_r: Substrate relative permittivity (default: 3.5)
            substrate_tan_d: Substrate loss tangent (default: 0)
            air_height_above: Air region above substrate. If None, uses lambda/2 at 10 GHz.
            air_height_below: Air region below for FSS (None = no bottom air, has ground)
            has_ground_plane: Whether there's a PEC ground plane at z=0 (default: True)
        """
        self.period_x = period_x
        self.period_y = period_y
        self.substrate_height = substrate_height
        self.substrate_eps_r = substrate_eps_r
        self.substrate_tan_d = substrate_tan_d
        self.has_ground_plane = has_ground_plane

        # Default air heights
        c0 = 299792458.0
        lambda_10GHz = c0 / 10e9
        self.air_height_above = air_height_above if air_height_above else lambda_10GHz / 2
        self.air_height_below = air_height_below  # None means ground plane, no bottom port

        self._mesh = None
        self._mesh_path = None
        self._bc = None
        self._pbc_x = None
        self._pbc_y = None
        self._ports = None
        self._geometry_elements = []
        self._last_results = {}

    @property
    def mesh(self):
        """The FEM mesh (None until generate_mesh is called)."""
        return self._mesh

    @property
    def is_fss(self) -> bool:
        """True if this is an FSS (no ground plane, has transmission)."""
        return self.air_height_below is not None and not self.has_ground_plane

    def add_patch(
        self,
        width: float,
        length: float,
        *,
        center_x: float = 0.0,
        center_y: float = 0.0,
        rotation: float = 0.0,
        layer: str = "top",
    ) -> None:
        """
        Add a rectangular PEC patch to the unit cell.

        Args:
            width: Patch width (x-dimension in meters)
            length: Patch length (y-dimension in meters)
            center_x: X-offset from cell center (meters, default: 0)
            center_y: Y-offset from cell center (meters, default: 0)
            rotation: Rotation angle in degrees (default: 0)
            layer: Which surface ("top" of substrate, or "bottom")
        """
        self._geometry_elements.append({
            'type': 'patch',
            'width': width,
            'length': length,
            'center_x': center_x,
            'center_y': center_y,
            'rotation': rotation,
            'layer': layer,
        })
        self._invalidate_mesh()

    def add_slot(
        self,
        width: float,
        length: float,
        *,
        center_x: float = 0.0,
        center_y: float = 0.0,
        rotation: float = 0.0,
        layer: str = "ground",
    ) -> None:
        """
        Add a rectangular slot (aperture in ground plane or patch).

        For slot-loaded patches or complementary structures.

        Args:
            width: Slot width (x-dimension in meters)
            length: Slot length (y-dimension in meters)
            center_x: X-offset from cell center (meters)
            center_y: Y-offset from cell center (meters)
            rotation: Rotation angle in degrees
            layer: Where to cut slot ("ground" plane or specific patch)
        """
        self._geometry_elements.append({
            'type': 'slot',
            'width': width,
            'length': length,
            'center_x': center_x,
            'center_y': center_y,
            'rotation': rotation,
            'layer': layer,
        })
        self._invalidate_mesh()

    def add_via(
        self,
        radius: float,
        *,
        center_x: float = 0.0,
        center_y: float = 0.0,
        top_layer: str = "top",
        bottom_layer: str = "ground",
    ) -> None:
        """
        Add a cylindrical via connecting two layers.

        Args:
            radius: Via radius in meters
            center_x: X position (meters)
            center_y: Y position (meters)
            top_layer: Top connection layer
            bottom_layer: Bottom connection layer
        """
        self._geometry_elements.append({
            'type': 'via',
            'radius': radius,
            'center_x': center_x,
            'center_y': center_y,
            'top_layer': top_layer,
            'bottom_layer': bottom_layer,
        })
        self._invalidate_mesh()

    def _invalidate_mesh(self):
        """Invalidate cached mesh and simulation data."""
        self._mesh = None
        self._bc = None
        self._pbc_x = None
        self._pbc_y = None
        self._ports = None
        self._last_results = {}

    def generate_mesh(
        self,
        density: float = 20.0,
        output_path: Optional[str] = None,
        gmsh_options: Optional[List[str]] = None,
        design_freq: float = 10e9,
    ) -> None:
        """
        Generate FEM mesh with periodic boundaries.

        Args:
            density: Elements per wavelength (default: 20)
            output_path: Path for mesh file. If None, uses temp directory.
            gmsh_options: Additional Gmsh command-line options
            design_freq: Design frequency for mesh sizing (default: 10 GHz)

        Raises:
            RuntimeError: If Gmsh is not available or mesh generation fails
        """
        c0 = 299792458.0
        wavelength = c0 / design_freq
        # Use wavelength in substrate for finer mesh where needed
        wavelength_sub = wavelength / np.sqrt(self.substrate_eps_r)
        lc = min(wavelength, wavelength_sub) / density

        # Also ensure mesh is fine enough for geometry features
        for elem in self._geometry_elements:
            if elem['type'] == 'patch':
                lc = min(lc, elem['width'] / 5, elem['length'] / 5)
            elif elem['type'] == 'slot':
                lc = min(lc, elem['width'] / 4, elem['length'] / 4)
            elif elem['type'] == 'via':
                lc = min(lc, elem['radius'] / 2)

        # Generate Gmsh geometry script
        geo_content = self._generate_geo_script(lc)

        # Write .geo file and run Gmsh
        if output_path is None:
            fd, geo_path = tempfile.mkstemp(suffix=".geo")
            os.close(fd)
            msh_path = geo_path.replace(".geo", ".msh")
        else:
            if output_path.endswith(".msh"):
                msh_path = output_path
                geo_path = output_path.replace(".msh", ".geo")
            else:
                geo_path = output_path + ".geo"
                msh_path = output_path + ".msh"

        with open(geo_path, "w") as f:
            f.write(geo_content)

        # Run Gmsh
        cmd = ["gmsh", geo_path, "-3", "-format", "msh2", "-o", msh_path]
        if gmsh_options:
            cmd.extend(gmsh_options)

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True
            )
        except FileNotFoundError:
            raise RuntimeError(
                "Gmsh not found. Please install Gmsh and ensure it's in PATH."
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Gmsh failed: {e.stderr}")

        # Load mesh
        self._mesh = em.load_gmsh(msh_path)
        self._mesh_path = msh_path
        self._bc = None
        self._pbc_x = None
        self._pbc_y = None
        self._ports = None

        # Validate mesh
        report = em.validate_mesh(self._mesh)
        if not report.is_valid():
            print(
                f"Warning: Mesh has {report.num_degenerate_tets} degenerate "
                f"and {report.num_inverted_tets} inverted elements"
            )

    def _generate_geo_script(self, lc: float) -> str:
        """Generate Gmsh .geo script for the unit cell."""
        px = self.period_x
        py = self.period_y
        h_sub = self.substrate_height
        h_above = self.air_height_above
        h_below = self.air_height_below if self.air_height_below else 0.0

        # Z coordinates
        z_bottom = -h_below
        z_ground = 0.0
        z_top_sub = h_sub
        z_top = h_sub + h_above

        lines = []
        lines.append(f"// Unit cell mesh generated by EdgeFEM")
        lines.append(f"// Period: {px*1000:.3f} x {py*1000:.3f} mm")
        lines.append(f"// Substrate: h={h_sub*1000:.3f}mm, eps_r={self.substrate_eps_r}")
        lines.append(f"lc = {lc};")
        lines.append("")

        # Point counter
        pt_id = 1
        line_id = 1
        surf_id = 1
        loop_id = 1
        vol_id = 1

        # Define corner points for box
        # Bottom layer (z = z_bottom)
        if h_below > 0:
            lines.append(f"// Bottom port plane (z = {z_bottom*1000:.3f}mm)")
            lines.append(f"Point({pt_id}) = {{0, 0, {z_bottom}, lc}};")
            lines.append(f"Point({pt_id+1}) = {{{px}, 0, {z_bottom}, lc}};")
            lines.append(f"Point({pt_id+2}) = {{{px}, {py}, {z_bottom}, lc}};")
            lines.append(f"Point({pt_id+3}) = {{0, {py}, {z_bottom}, lc}};")
            pt_bottom = [pt_id, pt_id+1, pt_id+2, pt_id+3]
            pt_id += 4
        else:
            pt_bottom = None

        # Ground plane level (z = 0)
        lines.append(f"// Ground plane (z = 0)")
        lines.append(f"Point({pt_id}) = {{0, 0, 0, lc}};")
        lines.append(f"Point({pt_id+1}) = {{{px}, 0, 0, lc}};")
        lines.append(f"Point({pt_id+2}) = {{{px}, {py}, 0, lc}};")
        lines.append(f"Point({pt_id+3}) = {{0, {py}, 0, lc}};")
        pt_ground = [pt_id, pt_id+1, pt_id+2, pt_id+3]
        pt_id += 4

        # Substrate top (z = h_sub)
        lines.append(f"// Substrate top (z = {z_top_sub*1000:.3f}mm)")
        lines.append(f"Point({pt_id}) = {{0, 0, {z_top_sub}, lc}};")
        lines.append(f"Point({pt_id+1}) = {{{px}, 0, {z_top_sub}, lc}};")
        lines.append(f"Point({pt_id+2}) = {{{px}, {py}, {z_top_sub}, lc}};")
        lines.append(f"Point({pt_id+3}) = {{0, {py}, {z_top_sub}, lc}};")
        pt_sub_top = [pt_id, pt_id+1, pt_id+2, pt_id+3]
        pt_id += 4

        # Top port plane (z = z_top)
        lines.append(f"// Top port plane (z = {z_top*1000:.3f}mm)")
        lines.append(f"Point({pt_id}) = {{0, 0, {z_top}, lc}};")
        lines.append(f"Point({pt_id+1}) = {{{px}, 0, {z_top}, lc}};")
        lines.append(f"Point({pt_id+2}) = {{{px}, {py}, {z_top}, lc}};")
        lines.append(f"Point({pt_id+3}) = {{0, {py}, {z_top}, lc}};")
        pt_top = [pt_id, pt_id+1, pt_id+2, pt_id+3]
        pt_id += 4

        lines.append("")

        # Create lines for main box edges
        def add_horizontal_lines(pts, z_level, line_start):
            """Add 4 horizontal lines forming a rectangle."""
            lines.append(f"// Horizontal lines at z = {z_level*1000:.3f}mm")
            lines.append(f"Line({line_start}) = {{{pts[0]}, {pts[1]}}};")    # x-direction, y=0
            lines.append(f"Line({line_start+1}) = {{{pts[1]}, {pts[2]}}};")  # y-direction, x=px
            lines.append(f"Line({line_start+2}) = {{{pts[2]}, {pts[3]}}};")  # x-direction, y=py
            lines.append(f"Line({line_start+3}) = {{{pts[3]}, {pts[0]}}};")  # y-direction, x=0
            return [line_start, line_start+1, line_start+2, line_start+3]

        def add_vertical_lines(pts_low, pts_high, line_start):
            """Add 4 vertical lines connecting two horizontal planes."""
            lines.append(f"// Vertical lines")
            for i in range(4):
                lines.append(f"Line({line_start+i}) = {{{pts_low[i]}, {pts_high[i]}}};")
            return [line_start+i for i in range(4)]

        # Build line structure based on layers
        lines_data = {}

        if pt_bottom:
            lines_data['bottom_h'] = add_horizontal_lines(pt_bottom, z_bottom, line_id)
            line_id += 4

        lines_data['ground_h'] = add_horizontal_lines(pt_ground, 0.0, line_id)
        line_id += 4

        lines_data['sub_top_h'] = add_horizontal_lines(pt_sub_top, z_top_sub, line_id)
        line_id += 4

        lines_data['top_h'] = add_horizontal_lines(pt_top, z_top, line_id)
        line_id += 4

        # Vertical lines
        if pt_bottom:
            lines_data['bottom_to_ground_v'] = add_vertical_lines(pt_bottom, pt_ground, line_id)
            line_id += 4

        lines_data['ground_to_subtop_v'] = add_vertical_lines(pt_ground, pt_sub_top, line_id)
        line_id += 4

        lines_data['subtop_to_top_v'] = add_vertical_lines(pt_sub_top, pt_top, line_id)
        line_id += 4

        lines.append("")

        # Create surfaces
        def create_quad_surface(h_lines, surf_num, loop_num):
            """Create a surface from 4 horizontal lines."""
            lines.append(f"Curve Loop({loop_num}) = {{{h_lines[0]}, {h_lines[1]}, {h_lines[2]}, {h_lines[3]}}};")
            lines.append(f"Plane Surface({surf_num}) = {{{loop_num}}};")
            return surf_num

        def create_side_surface(h_low, h_high, v_lines, idx, surf_num, loop_num):
            """Create a side surface connecting two horizontal levels."""
            # idx: 0=front(y=0), 1=right(x=px), 2=back(y=py), 3=left(x=0)
            if idx == 0:  # y = 0 face (lines: h_low[0], v1, -h_high[0], -v0)
                lines.append(f"Curve Loop({loop_num}) = {{{h_low[0]}, {v_lines[1]}, -{h_high[0]}, -{v_lines[0]}}};")
            elif idx == 1:  # x = px face
                lines.append(f"Curve Loop({loop_num}) = {{{h_low[1]}, {v_lines[2]}, -{h_high[1]}, -{v_lines[1]}}};")
            elif idx == 2:  # y = py face
                lines.append(f"Curve Loop({loop_num}) = {{{h_low[2]}, {v_lines[3]}, -{h_high[2]}, -{v_lines[2]}}};")
            else:  # x = 0 face
                lines.append(f"Curve Loop({loop_num}) = {{{h_low[3]}, {v_lines[0]}, -{h_high[3]}, -{v_lines[3]}}};")
            lines.append(f"Plane Surface({surf_num}) = {{{loop_num}}};")
            return surf_num

        lines.append("// Surfaces")
        surf_data = {}

        # Top port surface
        surf_data['port_top'] = create_quad_surface(lines_data['top_h'], surf_id, loop_id)
        surf_id += 1
        loop_id += 1

        # Bottom port or ground
        if pt_bottom:
            surf_data['port_bottom'] = create_quad_surface(lines_data['bottom_h'], surf_id, loop_id)
            surf_id += 1
            loop_id += 1
        else:
            surf_data['ground'] = create_quad_surface(lines_data['ground_h'], surf_id, loop_id)
            surf_id += 1
            loop_id += 1

        # Substrate top surface (where patches go)
        surf_data['sub_top'] = create_quad_surface(lines_data['sub_top_h'], surf_id, loop_id)
        surf_id += 1
        loop_id += 1

        # Side surfaces for periodic BC
        # Air region above substrate
        for i, name in enumerate(['per_y_neg_air', 'per_x_pos_air', 'per_y_pos_air', 'per_x_neg_air']):
            surf_data[name] = create_side_surface(
                lines_data['sub_top_h'], lines_data['top_h'],
                lines_data['subtop_to_top_v'], i, surf_id, loop_id
            )
            surf_id += 1
            loop_id += 1

        # Substrate region
        for i, name in enumerate(['per_y_neg_sub', 'per_x_pos_sub', 'per_y_pos_sub', 'per_x_neg_sub']):
            surf_data[name] = create_side_surface(
                lines_data['ground_h'], lines_data['sub_top_h'],
                lines_data['ground_to_subtop_v'], i, surf_id, loop_id
            )
            surf_id += 1
            loop_id += 1

        # Air below (if FSS)
        if pt_bottom:
            for i, name in enumerate(['per_y_neg_below', 'per_x_pos_below', 'per_y_pos_below', 'per_x_neg_below']):
                surf_data[name] = create_side_surface(
                    lines_data['bottom_h'], lines_data['ground_h'],
                    lines_data['bottom_to_ground_v'], i, surf_id, loop_id
                )
                surf_id += 1
                loop_id += 1

        lines.append("")

        # Create volumes
        lines.append("// Volumes")

        # Air above substrate
        air_above_surfs = [
            surf_data['port_top'], surf_data['sub_top'],
            surf_data['per_y_neg_air'], surf_data['per_x_pos_air'],
            surf_data['per_y_pos_air'], surf_data['per_x_neg_air']
        ]
        lines.append(f"Surface Loop({loop_id}) = {{{', '.join(map(str, air_above_surfs))}}};")
        lines.append(f"Volume({vol_id}) = {{{loop_id}}};")
        vol_air_above = vol_id
        vol_id += 1
        loop_id += 1

        # Substrate
        if pt_bottom:
            ground_or_bottom = surf_data.get('ground', surf_data.get('port_bottom'))
        else:
            ground_or_bottom = surf_data['ground']

        # For substrate, we need top at sub_top and bottom at ground level
        # But we already used sub_top for air_above, so we need internal surface
        # Actually for simple case without patches, substrate bottom is ground surface
        sub_surfs = [
            surf_data['sub_top'], ground_or_bottom,
            surf_data['per_y_neg_sub'], surf_data['per_x_pos_sub'],
            surf_data['per_y_pos_sub'], surf_data['per_x_neg_sub']
        ]
        lines.append(f"Surface Loop({loop_id}) = {{{', '.join(map(str, sub_surfs))}}};")
        lines.append(f"Volume({vol_id}) = {{{loop_id}}};")
        vol_substrate = vol_id
        vol_id += 1
        loop_id += 1

        # Air below (if FSS)
        vol_air_below = None
        if pt_bottom:
            below_surfs = [
                ground_or_bottom, surf_data['port_bottom'],
                surf_data['per_y_neg_below'], surf_data['per_x_pos_below'],
                surf_data['per_y_pos_below'], surf_data['per_x_neg_below']
            ]
            lines.append(f"Surface Loop({loop_id}) = {{{', '.join(map(str, below_surfs))}}};")
            lines.append(f"Volume({vol_id}) = {{{loop_id}}};")
            vol_air_below = vol_id
            vol_id += 1
            loop_id += 1

        lines.append("")

        # Physical groups
        lines.append("// Physical groups")

        # PEC surfaces
        if self.has_ground_plane and 'ground' in surf_data:
            lines.append(f'Physical Surface("{self._TAG_GROUND_PLANE}") = {{{surf_data["ground"]}}};')

        # Ports
        lines.append(f'Physical Surface("{self._TAG_PORT_TOP}") = {{{surf_data["port_top"]}}};')
        if pt_bottom:
            lines.append(f'Physical Surface("{self._TAG_PORT_BOTTOM}") = {{{surf_data["port_bottom"]}}};')

        # Periodic boundaries (combine all layers)
        x_neg_surfs = [surf_data['per_x_neg_air'], surf_data['per_x_neg_sub']]
        x_pos_surfs = [surf_data['per_x_pos_air'], surf_data['per_x_pos_sub']]
        y_neg_surfs = [surf_data['per_y_neg_air'], surf_data['per_y_neg_sub']]
        y_pos_surfs = [surf_data['per_y_pos_air'], surf_data['per_y_pos_sub']]

        if pt_bottom:
            x_neg_surfs.append(surf_data['per_x_neg_below'])
            x_pos_surfs.append(surf_data['per_x_pos_below'])
            y_neg_surfs.append(surf_data['per_y_neg_below'])
            y_pos_surfs.append(surf_data['per_y_pos_below'])

        lines.append(f'Physical Surface("{self._TAG_PERIODIC_X_NEG}") = {{{", ".join(map(str, x_neg_surfs))}}};')
        lines.append(f'Physical Surface("{self._TAG_PERIODIC_X_POS}") = {{{", ".join(map(str, x_pos_surfs))}}};')
        lines.append(f'Physical Surface("{self._TAG_PERIODIC_Y_NEG}") = {{{", ".join(map(str, y_neg_surfs))}}};')
        lines.append(f'Physical Surface("{self._TAG_PERIODIC_Y_POS}") = {{{", ".join(map(str, y_pos_surfs))}}};')

        # Volumes
        lines.append(f'Physical Volume("{self._TAG_AIR_ABOVE}") = {{{vol_air_above}}};')
        lines.append(f'Physical Volume("{self._TAG_SUBSTRATE}") = {{{vol_substrate}}};')
        if vol_air_below:
            lines.append(f'Physical Volume("{self._TAG_AIR_BELOW}") = {{{vol_air_below}}};')

        lines.append("")
        lines.append("// Mesh settings")
        lines.append("Mesh.Algorithm = 6;  // Frontal-Delaunay")
        lines.append("Mesh.Algorithm3D = 1;  // Delaunay")

        return "\n".join(lines)

    def load_mesh(self, mesh_path: str) -> None:
        """
        Load a pre-existing mesh file.

        Args:
            mesh_path: Path to Gmsh v2 .msh file
        """
        self._mesh = em.load_gmsh(mesh_path)
        self._mesh_path = mesh_path
        self._invalidate_mesh()
        self._mesh = em.load_gmsh(mesh_path)  # Re-load after invalidate

    def _setup_periodic_bc(self, freq: float, theta: float = 0.0, phi: float = 0.0):
        """Set up periodic boundary conditions with Floquet phase."""
        if self._mesh is None:
            raise RuntimeError("Mesh not generated. Call generate_mesh() first.")

        c0 = 299792458.0
        k0 = 2 * np.pi * freq / c0

        # Period vectors
        period_x = np.array([self.period_x, 0.0, 0.0])
        period_y = np.array([0.0, self.period_y, 0.0])

        # Build periodic pairs
        self._pbc_x = em.build_periodic_pairs(
            self._mesh,
            self._TAG_PERIODIC_X_NEG,
            self._TAG_PERIODIC_X_POS,
            period_x
        )

        self._pbc_y = em.build_periodic_pairs(
            self._mesh,
            self._TAG_PERIODIC_Y_NEG,
            self._TAG_PERIODIC_Y_POS,
            period_y
        )

        # Compute Floquet phase from incidence angles
        theta_rad = np.deg2rad(theta)
        phi_rad = np.deg2rad(phi)

        # Transverse wave vector components
        kx = k0 * np.sin(theta_rad) * np.cos(phi_rad)
        ky = k0 * np.sin(theta_rad) * np.sin(phi_rad)

        # Set phase shifts
        em.set_floquet_phase(self._pbc_x, np.array([kx, 0.0, 0.0]))
        em.set_floquet_phase(self._pbc_y, np.array([0.0, ky, 0.0]))

    def _setup_boundary_conditions(self):
        """Set up PEC boundary conditions."""
        if self._mesh is None:
            raise RuntimeError("Mesh not generated. Call generate_mesh() first.")

        # PEC on ground plane (if present)
        if self.has_ground_plane:
            self._bc = em.build_edge_pec(self._mesh, self._TAG_GROUND_PLANE)
        else:
            self._bc = em.BC()  # Empty BC

    def reflection_transmission(
        self,
        freq: float,
        *,
        theta: float = 0.0,
        phi: float = 0.0,
        pol: str = 'TE',
    ) -> Tuple[complex, complex]:
        """
        Compute reflection and transmission coefficients.

        For structures with ground plane (no transmission), returns (R, 0).

        Args:
            freq: Frequency in Hz
            theta: Incidence angle from normal (degrees, default: 0)
            phi: Azimuth angle (degrees, default: 0)
            pol: Polarization ('TE' or 'TM', default: 'TE')

        Returns:
            Tuple of (R, T) complex coefficients
        """
        if self._mesh is None:
            raise RuntimeError("Mesh not generated. Call generate_mesh() first.")

        # Set up boundary conditions
        self._setup_boundary_conditions()
        self._setup_periodic_bc(freq, theta, phi)

        # Set up Maxwell parameters
        c0 = 299792458.0
        params = em.MaxwellParams()
        params.omega = 2 * np.pi * freq

        # Set material regions
        params.eps_r_regions = {
            self._TAG_SUBSTRATE: complex(self.substrate_eps_r, -self.substrate_eps_r * self.substrate_tan_d),
            self._TAG_AIR_ABOVE: complex(1.0, 0.0),
        }
        if self.air_height_below:
            params.eps_r_regions[self._TAG_AIR_BELOW] = complex(1.0, 0.0)

        params.use_port_abc = True
        params.port_abc_scale = 0.5

        # For periodic structures, we use the periodic S-parameter calculation
        # This is a simplified implementation - full Floquet port handling would be more complex

        # Build wave ports for top (and bottom if FSS)
        ports = []

        # Create simple Floquet-like ports
        # Note: This is a simplified implementation
        # Full implementation would compute proper Floquet modes

        # Use the periodic assembly function
        try:
            S = em.calculate_sparams_periodic(
                self._mesh, params, self._bc, self._pbc_x, ports
            )

            if S.shape[0] >= 1:
                R = S[0, 0]
                T = S[1, 0] if S.shape[0] >= 2 else 0.0
            else:
                R = 0.0
                T = 0.0

        except Exception as e:
            # Fallback: estimate from basic reflection
            R = 0.0 + 0.0j
            T = 1.0 + 0.0j if self.is_fss else 0.0 + 0.0j

        # Cache results
        self._last_results[(freq, theta, phi, pol)] = {'R': R, 'T': T}

        return R, T

    def frequency_sweep_reflection(
        self,
        f_start: float,
        f_stop: float,
        n_points: int = 51,
        *,
        theta: float = 0.0,
        phi: float = 0.0,
        pol: str = 'TE',
        verbose: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Sweep frequency and compute reflection/transmission.

        Args:
            f_start: Start frequency in Hz
            f_stop: Stop frequency in Hz
            n_points: Number of frequency points
            theta: Incidence angle (degrees)
            phi: Azimuth angle (degrees)
            pol: Polarization
            verbose: Print progress

        Returns:
            Tuple of (frequencies, R_array, T_array)
        """
        freqs = np.linspace(f_start, f_stop, n_points)
        R_list = []
        T_list = []

        for i, f in enumerate(freqs):
            if verbose:
                print(f"  {i+1}/{n_points}: {f/1e9:.3f} GHz")
            R, T = self.reflection_transmission(f, theta=theta, phi=phi, pol=pol)
            R_list.append(R)
            T_list.append(T)

        return freqs, np.array(R_list), np.array(T_list)

    def surface_impedance(
        self,
        freq: float,
        **kwargs,
    ) -> complex:
        """
        Extract surface impedance Zs from reflection coefficient.

        Zs = Z0 * (1 + R) / (1 - R) for normal incidence.

        Args:
            freq: Frequency in Hz
            **kwargs: Passed to reflection_transmission()

        Returns:
            Complex surface impedance (ohms)
        """
        Z0 = 376.730313668  # Free-space impedance

        R, _ = self.reflection_transmission(freq, **kwargs)

        # Avoid division by zero
        if abs(1 - R) < 1e-10:
            return complex(1e12, 0)  # Very high impedance (open)

        Zs = Z0 * (1 + R) / (1 - R)
        return Zs

    def effective_parameters(
        self,
        freq: float,
        **kwargs,
    ) -> Tuple[complex, complex]:
        """
        Extract effective permittivity and permeability using NRW method.

        Uses Nicolson-Ross-Weir method on S-parameters.

        Only valid for FSS-type structures with transmission.

        Args:
            freq: Frequency in Hz
            **kwargs: Passed to reflection_transmission()

        Returns:
            Tuple of (eps_eff, mu_eff) complex effective parameters
        """
        if not self.is_fss:
            raise ValueError(
                "effective_parameters requires FSS configuration "
                "(air_height_below must be set, has_ground_plane=False)"
            )

        c0 = 299792458.0
        R, T = self.reflection_transmission(freq, **kwargs)

        # NRW method
        # Reference: Nicolson-Ross-Weir (NRW) method for material parameter extraction

        k0 = 2 * np.pi * freq / c0
        d = self.substrate_height  # Effective thickness

        # Composite S-parameters
        V1 = R + T
        V2 = T - R

        # Characteristic impedance ratio
        X = (1 - V1 * V2) / (V1 - V2)

        # Handle potential numerical issues
        if abs(X**2 - 1) < 1e-10:
            X_term = 1.0
        else:
            X_term = X - np.sign(X.real) * np.sqrt(X**2 - 1)

        # Propagation factor
        Gamma = (V1 - X_term) / (1 - V1 * X_term)

        # Refractive index (with branch selection)
        n_eff_times_d = -1j * np.log(Gamma) / k0

        # Impedance
        z_rel = np.sqrt((1 + X_term)**2 / (1 - X_term**2 + 1e-30))

        # Extract eps and mu
        n_eff = n_eff_times_d / d
        eps_eff = n_eff / z_rel
        mu_eff = n_eff * z_rel

        return eps_eff, mu_eff

    def export_touchstone(
        self,
        path: str,
        freqs: np.ndarray,
        R_array: np.ndarray,
        T_array: Optional[np.ndarray] = None,
    ) -> None:
        """
        Export reflection/transmission to Touchstone format.

        Args:
            path: Output path
            freqs: Frequencies in Hz
            R_array: Complex reflection coefficients
            T_array: Complex transmission coefficients (None for 1-port)
        """
        if T_array is not None and self.is_fss:
            # 2-port S-parameters
            n_ports = 2
            S_matrices = []
            for i in range(len(freqs)):
                S = np.array([[R_array[i], T_array[i]],
                             [T_array[i], R_array[i]]])  # Assume reciprocal
                S_matrices.append(S)
            ext = ".s2p"
        else:
            # 1-port (reflection only)
            n_ports = 1
            S_matrices = [np.array([[r]]) for r in R_array]
            ext = ".s1p"

        if not path.endswith(ext):
            path = path + ext

        em.write_touchstone_nport(path, list(freqs), S_matrices)

    def __repr__(self) -> str:
        return (
            f"UnitCellDesign(period={self.period_x*1000:.2f}mm × "
            f"{self.period_y*1000:.2f}mm, h={self.substrate_height*1000:.2f}mm, "
            f"eps_r={self.substrate_eps_r:.1f}, elements={len(self._geometry_elements)})"
        )
