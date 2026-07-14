#pragma once

#include <unordered_set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "edgefem/mesh.hpp"
#include "edgefem/ports/port_eigensolve.hpp"

namespace edgefem {

struct PortSurfaceMesh {
  Mesh mesh;
  std::vector<int> volume_tri_indices; // index into volume mesh tris
};

PortSurfaceMesh extract_surface_mesh(const Mesh &volume_mesh, int surface_tag);

struct WavePort {
  int surface_tag = 0;
  PortMode mode;
  std::vector<int> edges;   // global edge indices participating in the port
  Eigen::VectorXcd weights; // line-integral weights for each edge
};

WavePort build_wave_port(const Mesh &volume_mesh,
                         const PortSurfaceMesh &surface, const PortMode &mode);

/// Populate TE10 mode field on a surface mesh using analytical formula.
/// The surface mesh must be a rectangular cross-section in the XY plane.
/// @param surface The port surface mesh (extracted from volume mesh)
/// @param port Rectangular waveguide port dimensions
/// @param mode PortMode with fc, omega, etc. populated (field will be filled)
void populate_te10_field(const PortSurfaceMesh &surface,
                         const RectWaveguidePort &port, PortMode &mode);

/// Build wave port using 3D FEM eigenvector as weights.
/// This provides better coupling to the actual FEM mode than analytical
/// weights.
/// @param volume_mesh The full 3D mesh
/// @param surface The port surface mesh
/// @param eigenvector Full 3D FEM eigenvector (size = num_edges)
/// @param mode PortMode with Z0, beta, etc. populated
/// @param bc Boundary conditions (to identify PEC edges)
/// @return WavePort with eigenvector-based weights
WavePort build_wave_port_from_eigenvector(
    const Mesh &volume_mesh, const PortSurfaceMesh &surface,
    const Eigen::VectorXd &eigenvector, const PortMode &mode,
    const std::unordered_set<int> &pec_edges);

/// Compute 3D FEM eigenvector for a TE mode with given cutoff wavenumber.
/// Solves generalized eigenvalue problem (K, M) to find mode closest to
/// target_kc_sq.
/// @param mesh The full 3D mesh
/// @param pec_edges Set of PEC edge indices (Dirichlet BC)
/// @param target_kc_sq Target cutoff wavenumber squared (e.g., (π/a)² for TE10)
/// @return Eigenvector (size = num_edges, zero for PEC edges)
Eigen::VectorXd compute_te_eigenvector(const Mesh &mesh,
                                       const std::unordered_set<int> &pec_edges,
                                       double target_kc_sq);

/// Solve the 2D discrete TE port mode directly on the port face, in the
/// same tangential edge basis the 3D system uses there.
/// Solves the generalized eigenproblem K_s e = kc² M_s e restricted to the
/// port's free (non-PEC) edges, where K_s is the surface curl-curl matrix
/// and M_s the surface mass matrix, and returns the eigenvector whose kc²
/// is closest to target_kc_sq (gradient-nullspace modes kc²≈0 are skipped).
/// The discrete kc² is written to kc_sq_out — use it (not the analytical
/// value) to compute β for the port ABC so profile and impedance are
/// mutually consistent.
/// Dense solve; intended for port cross-sections (≲ a few thousand edges).
/// @return Eigenvector in GLOBAL edge indexing (zero off the port)
Eigen::VectorXd solve_port_mode_2d(const Mesh &mesh, int surface_tag,
                                   const std::unordered_set<int> &pec_edges,
                                   double target_kc_sq, double &kc_sq_out);

/// Build a wave port from the 2D discrete port-face eigenmode (see
/// solve_port_mode_2d). Replaces the 3D-cavity-eigenvector approach: the
/// weights are the discrete port mode itself, and mode.kc is set to the
/// DISCRETE cutoff wavenumber so calculate_sparams_eigenmode uses a
/// consistent β. Much cheaper than compute_te_eigenvector (2D dense solve
/// on port edges instead of a 3D dense eigensolve).
WavePort build_wave_port_2d(const Mesh &mesh, int surface_tag,
                            const PortMode &mode,
                            const std::unordered_set<int> &pec_edges,
                            double target_kc_sq);

/// Assemble the port surface mass matrix
///   M_s[i][j] = ∫_S (n̂×N_i)·(n̂×N_j) dS = ∫_S N_i·N_j dS
/// over the triangles of the port surface, in GLOBAL edge indexing
/// (size num_edges × num_edges, nonzeros only on port edges). Rows/columns
/// of Dirichlet (PEC) edges are omitted. This is the operator required by
/// the first-order modal ABC  n̂×(∇×E) + jβ E_t = 2jβ E_inc,t  in weak form.
/// @param mesh The full 3D mesh
/// @param surface_tag Physical tag of the port surface
/// @param dirichlet_edges PEC edge indices to exclude
Eigen::SparseMatrix<double>
assemble_port_surface_mass(const Mesh &mesh, int surface_tag,
                           const std::unordered_set<int> &dirichlet_edges);

} // namespace edgefem
