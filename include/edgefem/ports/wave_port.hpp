#pragma once

#include <unordered_set>
#include <vector>

#include <Eigen/Core>

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
  std::vector<int> edges;      // global edge indices participating in the port
  Eigen::VectorXcd weights;    // line-integral weights for each edge
};

WavePort build_wave_port(const Mesh &volume_mesh, const PortSurfaceMesh &surface,
                         const PortMode &mode);

/// Populate TE10 mode field on a surface mesh using analytical formula.
/// The surface mesh must be a rectangular cross-section in the XY plane.
/// @param surface The port surface mesh (extracted from volume mesh)
/// @param port Rectangular waveguide port dimensions
/// @param mode PortMode with fc, omega, etc. populated (field will be filled)
void populate_te10_field(const PortSurfaceMesh &surface,
                         const RectWaveguidePort &port,
                         PortMode &mode);

/// Build wave port using 3D FEM eigenvector as weights.
/// This provides better coupling to the actual FEM mode than analytical weights.
/// @param volume_mesh The full 3D mesh
/// @param surface The port surface mesh
/// @param eigenvector Full 3D FEM eigenvector (size = num_edges)
/// @param mode PortMode with Z0, beta, etc. populated
/// @param bc Boundary conditions (to identify PEC edges)
/// @return WavePort with eigenvector-based weights
WavePort build_wave_port_from_eigenvector(const Mesh &volume_mesh,
                                           const PortSurfaceMesh &surface,
                                           const Eigen::VectorXd &eigenvector,
                                           const PortMode &mode,
                                           const std::unordered_set<int> &pec_edges);

/// Compute 3D FEM eigenvector for a TE mode with given cutoff wavenumber.
/// Solves generalized eigenvalue problem (K, M) to find mode closest to target_kc_sq.
/// @param mesh The full 3D mesh
/// @param pec_edges Set of PEC edge indices (Dirichlet BC)
/// @param target_kc_sq Target cutoff wavenumber squared (e.g., (π/a)² for TE10)
/// @return Eigenvector (size = num_edges, zero for PEC edges)
Eigen::VectorXd compute_te_eigenvector(const Mesh &mesh,
                                        const std::unordered_set<int> &pec_edges,
                                        double target_kc_sq);

} // namespace edgefem

