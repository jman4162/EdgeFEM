import numpy as np
from vectorem_cad.sdk.project import Project

def run_sdk_demo():
    """
    Demonstrates the use of the high-level Python SDK to run a simulation.
    """
    # 1. Create a project
    proj = Project(name="waveguide_sdk_demo")

    # 2. Load the mesh
    proj.load_mesh_from_gmsh("examples/cube_cavity.msh")

    # 3. Define boundaries and ports
    proj.set_pec_boundary(tag=1)
    proj.add_lumped_port(edge_id=10, z0=50.0)
    proj.add_lumped_port(edge_id=20, z0=50.0)

    # 4. Define the frequency sweep
    freqs = np.linspace(8e9, 12e9, 11)

    # 5. Run the simulation
    print("\\nRunning frequency sweep via SDK...")
    for f in freqs:
        results = proj.solve_s_params(f)
        print(f"f={f/1e9:.2f} GHz, S11={results.s11:.4f}")

if __name__ == "__main__":
    run_sdk_demo()
