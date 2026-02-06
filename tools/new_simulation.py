#!/usr/bin/env python3
"""
EdgeFEM New Simulation Setup

Creates a new simulation directory with geometry, mesh, and run script.

Usage:
    python tools/new_simulation.py waveguide my_project
    python tools/new_simulation.py patch my_antenna --freq 2.45e9
    python tools/new_simulation.py unit_cell my_array --cell-size 15
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from tools.geometry import RectWaveguide, PatchAntenna, UnitCell, DipoleAntenna


def create_run_script(sim_dir: Path, geo_file: str, sim_type: str):
    """Create a Python run script for the simulation."""
    script = f'''#!/usr/bin/env python3
"""
EdgeFEM Simulation: {sim_dir.name}
Type: {sim_type}

Run with: python run_simulation.py
Can also be run from any directory.
"""

import sys
from pathlib import Path
import numpy as np

# Get the directory where this script lives (simulation directory)
SCRIPT_DIR = Path(__file__).parent.resolve()

# Add EdgeFEM to path
edgefem_root = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(edgefem_root / "build" / "python"))

import pyedgefem as em


def main():
    # Configuration - use absolute paths based on script location
    geo_file = SCRIPT_DIR / "{geo_file}"
    msh_file = geo_file.with_suffix(".msh")

    # Mesh the geometry (if not already done)
    if not msh_file.exists():
        import subprocess
        print(f"Meshing {{geo_file.name}}...")
        # Use -format msh2 for Gmsh v2 format (required by EdgeFEM)
        result = subprocess.run(["gmsh", str(geo_file), "-3", "-format", "msh2", "-o", str(msh_file)],
                                capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Meshing failed: {{result.stderr}}")
            return

    # Load mesh
    print(f"Loading mesh from {{msh_file.name}}...")
    mesh = em.load_gmsh(str(msh_file))
    print(f"  Nodes: {{mesh.num_nodes()}}")
    print(f"  Tets: {{mesh.num_tets()}}")
    print(f"  Tris: {{mesh.num_tris()}}")

    # Validate mesh
    quality = em.validate_mesh(mesh)
    if quality.is_valid():
        print("  Mesh quality: OK")
    else:
        print(f"  Warning: {{quality.num_degenerate_tets}} degenerate, "
              f"{{quality.num_inverted_tets}} inverted tets")

    # Set up boundary conditions
    # Adjust PEC tag based on your geometry
    bc = em.build_edge_pec(mesh, 1)
    print(f"  Applied PEC boundary")

    # Frequency setup
    freq = 10e9  # Modify as needed
    omega = 2 * np.pi * freq

    params = em.MaxwellParams()
    params.omega = omega
    params.eps_r = 1.0 + 0j
    params.mu_r = 1.0 + 0j

    # Solve
    print(f"\\nSolving at f = {{freq/1e9:.2f}} GHz...")

    # For simulations without ports, just assemble and solve
    asm = em.assemble_maxwell(mesh, params, bc, [], -1)

    opts = em.SolveOptions()
    opts.tolerance = 1e-10
    opts.max_iterations = 10000

    result = em.solve_linear(asm.A, asm.b, opts)

    print(f"  Converged: {{result.converged}}")
    print(f"  Iterations: {{result.iters}}")
    print(f"  Residual: {{result.residual:.2e}}")

    # TODO: Add your post-processing here
    # - S-parameter extraction
    # - Far-field patterns
    # - Field visualization

    print("\\nSimulation complete!")


if __name__ == "__main__":
    main()
'''
    (sim_dir / "run_simulation.py").write_text(script)


def create_waveguide_simulation(sim_dir: Path, args):
    """Create waveguide simulation."""
    wg = RectWaveguide(
        a=args.width / 1000 if args.width else 0.02286,
        b=args.height / 1000 if args.height else 0.01016,
        length=args.length / 1000 if args.length else 0.05
    )
    geo_file = "waveguide.geo"
    wg.generate(str(sim_dir / geo_file))
    create_run_script(sim_dir, geo_file, "waveguide")
    return geo_file


def create_patch_simulation(sim_dir: Path, args):
    """Create patch antenna simulation."""
    # Calculate patch dimensions from frequency if provided
    c0 = 299792458
    eps_r = 4.4

    if args.freq:
        wavelength = c0 / args.freq
        effective_wavelength = wavelength / np.sqrt(eps_r)
        patch_length = effective_wavelength / 2 * 0.9  # Approximate
        patch_width = patch_length * 1.3
    else:
        patch_width = 0.038
        patch_length = 0.029

    patch = PatchAntenna(
        patch_width=patch_width,
        patch_length=patch_length,
        substrate_height=args.substrate_h / 1000 if args.substrate_h else 0.00157,
        substrate_eps_r=eps_r
    )
    geo_file = "patch_antenna.geo"
    patch.generate(str(sim_dir / geo_file))
    create_run_script(sim_dir, geo_file, "patch_antenna")
    return geo_file


def create_unit_cell_simulation(sim_dir: Path, args):
    """Create unit cell simulation."""
    cell_size = args.cell_size / 1000 if args.cell_size else 0.015

    unit = UnitCell(
        cell_size_x=cell_size,
        cell_size_y=cell_size,
        element_width=cell_size * 0.5,
        element_length=cell_size * 0.65,
        substrate_height=args.substrate_h / 1000 if args.substrate_h else 0.00157
    )
    geo_file = "unit_cell.geo"
    unit.generate(str(sim_dir / geo_file))
    create_run_script(sim_dir, geo_file, "unit_cell")
    return geo_file


def create_dipole_simulation(sim_dir: Path, args):
    """Create dipole antenna simulation."""
    freq = args.freq if args.freq else 1e9
    c0 = 299792458
    wavelength = c0 / freq

    dipole = DipoleAntenna(
        length=wavelength / 2,
        frequency=freq
    )
    geo_file = "dipole.geo"
    dipole.generate(str(sim_dir / geo_file))
    create_run_script(sim_dir, geo_file, "dipole")
    return geo_file


def main():
    parser = argparse.ArgumentParser(
        description="Create a new EdgeFEM simulation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python tools/new_simulation.py waveguide my_wg_project
    python tools/new_simulation.py patch my_patch --freq 2.45e9
    python tools/new_simulation.py unit_cell my_array --cell-size 15
    python tools/new_simulation.py dipole my_dipole --freq 300e6

Simulation types:
    waveguide   - Rectangular waveguide with wave ports
    patch       - Microstrip patch antenna on substrate
    unit_cell   - Periodic unit cell for phased array
    dipole      - Half-wave dipole in free space
"""
    )

    parser.add_argument("type", choices=["waveguide", "patch", "unit_cell", "dipole"],
                        help="Type of simulation to create")
    parser.add_argument("name", help="Name for the simulation directory")
    parser.add_argument("--freq", type=float, help="Design frequency in Hz")
    parser.add_argument("--width", type=float, help="Width in mm (waveguide)")
    parser.add_argument("--height", type=float, help="Height in mm (waveguide)")
    parser.add_argument("--length", type=float, help="Length in mm")
    parser.add_argument("--cell-size", type=float, help="Cell size in mm (unit_cell)")
    parser.add_argument("--substrate-h", type=float, help="Substrate height in mm")
    parser.add_argument("--mesh", action="store_true", help="Run gmsh to create mesh")

    args = parser.parse_args()

    # Create simulation directory
    sim_dir = project_root / "simulations" / args.name
    sim_dir.mkdir(parents=True, exist_ok=True)
    print(f"Creating simulation in: {sim_dir}")

    # Generate geometry based on type
    if args.type == "waveguide":
        geo_file = create_waveguide_simulation(sim_dir, args)
    elif args.type == "patch":
        geo_file = create_patch_simulation(sim_dir, args)
    elif args.type == "unit_cell":
        geo_file = create_unit_cell_simulation(sim_dir, args)
    elif args.type == "dipole":
        geo_file = create_dipole_simulation(sim_dir, args)

    # Optionally mesh
    if args.mesh:
        geo_path = sim_dir / geo_file
        msh_path = geo_path.with_suffix(".msh")
        print(f"Meshing {geo_file}...")
        # Use -format msh2 for Gmsh v2 format (required by EdgeFEM)
        result = subprocess.run(
            ["gmsh", str(geo_path), "-3", "-format", "msh2", "-o", str(msh_path)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"Created: {msh_path}")
        else:
            print(f"Meshing failed: {result.stderr}")

    print()
    print("Next steps:")
    print(f"  1. cd simulations/{args.name}")
    if not args.mesh:
        print(f"  2. gmsh {geo_file} -3 -format msh2 -o {geo_file.replace('.geo', '.msh')}")
        print(f"  3. python run_simulation.py")
    else:
        print(f"  2. python run_simulation.py")


if __name__ == "__main__":
    main()
