#pragma once

#include <complex>
#include <vector>

#include <Eigen/Core>

#include "edgefem/mesh.hpp"

namespace edgefem {

using VecC = Eigen::VectorXcd;

struct HuygensSurfaceData {
  std::vector<Eigen::Vector3d> r;       // Triangle centroids
  std::vector<Eigen::Vector3d> n;       // Outward normals
  std::vector<Eigen::Vector3cd> E_tan;  // Tangential E
  std::vector<Eigen::Vector3cd> H_tan;  // Tangential H
  std::vector<double> area;             // Triangle areas
};

/// Extract Huygens surface data from a FEM solution on a tagged surface.
/// Computes tangential E and H at each triangle centroid by finding the parent
/// tetrahedron containing each surface triangle and evaluating the Whitney
/// basis functions.
///
/// The output is directly compatible with stratton_chu_2d() and
/// stratton_chu_3d() for near-to-far field transformation.
///
/// @param mesh Volume mesh with surface triangles
/// @param solution Edge-element solution vector (complex)
/// @param surface_tag Physical tag of the Huygens surface
/// @param omega Angular frequency (rad/s)
/// @param mu_r Relative permeability (default: 1.0)
/// @return HuygensSurfaceData with tangential fields at triangle centroids
HuygensSurfaceData extract_huygens_surface(const Mesh &mesh,
                                            const VecC &solution,
                                            int surface_tag, double omega,
                                            std::complex<double> mu_r = 1.0);

} // namespace edgefem
