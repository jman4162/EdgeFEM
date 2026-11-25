#pragma once

#include <complex>
#include <optional>
#include <vector>

#include <Eigen/Core>

#include "edgefem/mesh.hpp"

namespace edgefem {

using VecC = Eigen::VectorXcd;

struct HuygensSurfaceData {
  std::vector<Eigen::Vector3d> r;      // Triangle centroids
  std::vector<Eigen::Vector3d> n;      // Outward normals
  std::vector<Eigen::Vector3cd> E_tan; // Tangential E
  std::vector<Eigen::Vector3cd> H_tan; // Tangential H
  std::vector<double> area;            // Triangle areas
};

/// Extract Huygens surface data from a FEM solution on a tagged surface.
/// For each triangle with `surface_tag`, E is evaluated at the centroid from
/// the Whitney interpolant and H = -curl E / (j omega mu0 mu_r) (time
/// convention exp(+j omega t)), each averaged over the one or two tets that
/// share the triangle, then projected onto the triangle plane.
///
/// Normals: with `interior_point`, each normal points away from it (use a
/// point inside the closed Huygens surface, or on the source side of an open
/// one). Without it, a boundary triangle's normal points out of its one
/// parent tet; on an interior surface the side is arbitrary, so pass
/// `interior_point` there.
///
/// The output feeds stratton_chu_2d() and stratton_chu_3d() directly.
///
/// @param mesh Volume mesh with surface triangles
/// @param solution Edge-element solution vector (complex)
/// @param surface_tag Physical tag of the Huygens surface
/// @param omega Angular frequency (rad/s)
/// @param mu_r Relative permeability at the surface (default: 1.0)
/// @param interior_point Reference point the normals point away from
/// @return HuygensSurfaceData with tangential fields at triangle centroids
HuygensSurfaceData
extract_huygens_surface(const Mesh &mesh, const VecC &solution, int surface_tag,
                        double omega, std::complex<double> mu_r = 1.0,
                        std::optional<Eigen::Vector3d> interior_point = {});

} // namespace edgefem
