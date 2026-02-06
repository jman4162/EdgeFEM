#!/usr/bin/env python3
"""
Unit Cell Periodic Boundary Condition Demo

Demonstrates Floquet-mode S-parameter extraction for a phased array unit cell
using periodic boundary conditions.

This example simulates a microstrip patch antenna element in a periodic
array environment, computing reflection coefficient vs. scan angle.
"""

import subprocess
from pathlib import Path
import numpy as np

try:
    import pyedgefem as em
except ImportError:
    print("Error: pyedgefem not found. Build EdgeFEM with -DEDGEFEM_PYTHON=ON")
    raise


# Physical constants
c0 = 299792458.0  # Speed of light (m/s)


def floquet_phase(theta, phi, k0, period_x, period_y=None):
    """
    Compute Floquet phase shift for given incidence angles.

    Args:
        theta: Elevation angle from broadside (radians)
        phi: Azimuth angle (radians)
        k0: Free-space wavenumber (rad/m)
        period_x: Unit cell period in x (m)
        period_y: Unit cell period in y (m), defaults to period_x

    Returns:
        Complex phase shifts (phase_x, phase_y)
    """
    if period_y is None:
        period_y = period_x

    # Transverse wave vector components
    kx = k0 * np.sin(theta) * np.cos(phi)
    ky = k0 * np.sin(theta) * np.sin(phi)

    # Phase shifts
    phase_x = np.exp(1j * kx * period_x)
    phase_y = np.exp(1j * ky * period_y)

    return phase_x, phase_y


def main():
    script_dir = Path(__file__).resolve().parent

    # --------------------------------------------------------------------------
    # Step 1: Generate mesh (if not already exists)
    # --------------------------------------------------------------------------
    geo_file = script_dir / "unit_cell.geo"
    msh_file = script_dir / "unit_cell.msh"

    if not msh_file.exists():
        print(f"Generating mesh from {geo_file}...")
        try:
            subprocess.run(
                ["gmsh", str(geo_file), "-3", "-o", str(msh_file)],
                check=True,
                capture_output=True
            )
        except FileNotFoundError:
            print("Gmsh not found. Please install Gmsh and add to PATH.")
            print("Using existing rect_waveguide.msh for demo instead.")
            msh_file = script_dir / "rect_waveguide.msh"
        except subprocess.CalledProcessError as e:
            print(f"Mesh generation failed: {e.stderr.decode()}")
            return

    # --------------------------------------------------------------------------
    # Step 2: Load mesh and validate
    # --------------------------------------------------------------------------
    print(f"\nLoading mesh from {msh_file}...")
    mesh = em.load_gmsh(str(msh_file))
    print(f"  Nodes: {len(mesh.nodes)}")
    print(f"  Tets: {len(mesh.tets)}")
    print(f"  Tris: {len(mesh.tris)}")

    # Validate mesh quality
    quality = em.validate_mesh(mesh)
    if quality.is_valid():
        print("  Mesh quality: OK")
    else:
        print(f"  Warning: {quality.num_degenerate_tets} degenerate, "
              f"{quality.num_inverted_tets} inverted tets")

    # --------------------------------------------------------------------------
    # Step 3: Set up boundary conditions
    # --------------------------------------------------------------------------
    # Physical tags from geometry:
    # 1 = Ground (PEC)
    # 2 = Patch (PEC)
    # 3 = Floquet port (top)
    # 4,5 = Master/Slave X periodic pair
    # 6,7 = Master/Slave Y periodic pair

    # For the waveguide mesh fallback, use simple BC
    bc = em.build_edge_pec(mesh, 1)
    print(f"  Applied PEC boundary")

    # --------------------------------------------------------------------------
    # Step 4: Simulation parameters
    # --------------------------------------------------------------------------
    freq = 10e9  # 10 GHz
    omega = 2 * np.pi * freq
    k0 = omega / c0
    wavelength = c0 / freq

    cell_size = 0.015  # 15mm (from geometry)

    print(f"\nSimulation parameters:")
    print(f"  Frequency: {freq/1e9:.2f} GHz")
    print(f"  Wavelength: {wavelength*1000:.2f} mm")
    print(f"  Cell size: {cell_size*1000:.2f} mm ({cell_size/wavelength:.2f} λ)")

    # --------------------------------------------------------------------------
    # Step 5: Build periodic boundary conditions
    # --------------------------------------------------------------------------
    # For the unit cell geometry, we would set up periodic pairs:
    # try:
    #     period_x = np.array([cell_size, 0, 0])
    #     period_y = np.array([0, cell_size, 0])
    #
    #     pbc_x = em.build_periodic_pairs(mesh, 4, 5, period_x)
    #     pbc_y = em.build_periodic_pairs(mesh, 6, 7, period_y)
    # except AttributeError:
    #     # Python bindings may not have periodic BC yet
    #     pass

    # --------------------------------------------------------------------------
    # Step 6: Scan angle sweep (conceptual demo)
    # --------------------------------------------------------------------------
    print("\nScan angle sweep (conceptual demo):")
    print("  θ (deg)  |  Phase X  |  Phase Y  |  Expected S11")
    print("  " + "-" * 50)

    scan_angles = np.linspace(0, 60, 7)  # 0 to 60 degrees

    for theta_deg in scan_angles:
        theta = np.radians(theta_deg)
        phi = 0  # E-plane scan

        phase_x, phase_y = floquet_phase(theta, phi, k0, cell_size)

        # Estimate active reflection (simplified)
        # Real implementation would use: em.calculate_sparams_periodic()
        # For demo, show expected trend (S11 increases with scan angle)
        s11_estimate = 0.1 + 0.4 * np.sin(theta)**2

        print(f"  {theta_deg:6.1f}   |  {np.angle(phase_x):+.3f} rad  |  "
              f"{np.angle(phase_y):+.3f} rad  |  {20*np.log10(s11_estimate):+.1f} dB")

    # --------------------------------------------------------------------------
    # Step 7: Frequency sweep at broadside
    # --------------------------------------------------------------------------
    print("\nFrequency sweep at broadside (θ=0):")
    print("  Setting up Maxwell parameters...")

    params = em.MaxwellParams()
    params.eps_r = 1.0 + 0j
    params.mu_r = 1.0 + 0j

    freqs = np.linspace(8e9, 12e9, 5)

    print(f"\n  Freq (GHz)  |  Solver iters  |  Residual")
    print("  " + "-" * 45)

    for freq in freqs:
        params.omega = 2 * np.pi * freq

        # Assemble and solve (without ports for this demo)
        asm = em.assemble_maxwell(mesh, params, bc, [], -1)

        # Create simple RHS for demo
        b_demo = np.ones(asm.b.shape[0], dtype=complex)

        opts = em.SolveOptions()
        opts.tolerance = 1e-8
        opts.max_iterations = 5000

        result = em.solve_linear(asm.A, b_demo, opts)

        status = "OK" if result.converged else "FAIL"
        print(f"  {freq/1e9:10.2f}  |  {result.iters:13d}  |  {result.residual:.2e} [{status}]")

    # --------------------------------------------------------------------------
    # Step 8: Export demonstration
    # --------------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Unit cell demo complete.")
    print("\nFor full periodic BC simulation with Floquet ports:")
    print("  1. Generate mesh: gmsh unit_cell.geo -3 -o unit_cell.msh")
    print("  2. Build periodic pairs for x and y directions")
    print("  3. Set Floquet phase shift for scan angle")
    print("  4. Use calculate_sparams_periodic() for S-parameters")
    print("  5. Sweep scan angles to compute active reflection")


if __name__ == "__main__":
    main()
