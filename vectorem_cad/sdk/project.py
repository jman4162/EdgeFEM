import numpy as np
import pyvectorem as em

class Project:
    """
    High-level class for managing a VectorEM simulation project.
    """
    def __init__(self, name="default"):
        self.name = name
        self.mesh = None
        self.bcs = []
        self.ports = []
        self.params = em.MaxwellParams()

    def load_mesh_from_gmsh(self, filepath):
        """Loads a mesh from a Gmsh .msh file."""
        print(f"Loading mesh from {filepath}...")
        self.mesh = em.load_gmsh(filepath)
        print("Mesh loaded successfully.")

    def set_pec_boundary(self, tag=1):
        """Defines a PEC boundary condition from a physical tag."""
        if not self.mesh:
            raise RuntimeError("Mesh must be loaded before setting BCs.")
        bc = em.build_edge_pec(self.mesh, tag)
        self.bcs.append(bc) # This is not quite right, BCs should be merged
        print(f"PEC boundary defined for tag {tag}.")

    def add_lumped_port(self, edge_id, z0=50.0):
        """Adds a simplified edge-based wave port to the simulation."""
        port = em.WavePort()
        port.edges = [edge_id]
        port.weights = np.array([1.0 + 0.0j])
        mode = em.PortMode()
        mode.Z0 = z0
        port.mode = mode
        self.ports.append(port)
        print(f"Added edge port on edge {edge_id} with Z0={z0} Ohms.")

    def solve_s_params(self, freq):
        """Solves for the S-parameters at a given frequency."""
        if not self.mesh:
            raise RuntimeError("Mesh must be loaded before solving.")

        self.params.omega = 2 * 3.1415926535 * freq

        # This assumes a single BC object for all PECs.
        # A more robust implementation would merge BCs.
        bc = self.bcs[0] if self.bcs else em.BC()

        print(f"Solving at {freq/1e9:.2f} GHz...")
        S = em.calculate_sparams(self.mesh, self.params, bc, self.ports)
        print("Solve complete.")

        return Results(S)


class Results:
    """
    A container for simulation results.
    """
    def __init__(self, s_matrix):
        self.s_matrix = s_matrix

    def __repr__(self):
        return f"Results(S-matrix shape: {self.s_matrix.shape})"

    @property
    def s11(self):
        return self.s_matrix[0, 0]

    @property
    def s21(self):
        return self.s_matrix[1, 0]

    @property
    def s12(self):
        return self.s_matrix[0, 1]

    @property
    def s22(self):
        return self.s_matrix[1, 1]
