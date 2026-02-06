#pragma once

#include <vector>

#include "edgefem/mesh.hpp"

namespace edgefem {

struct MeshQualityReport {
  int num_degenerate_tets = 0;
  int num_inverted_tets = 0;
  double min_volume = 0.0;
  double max_aspect_ratio = 0.0;
  std::vector<int> bad_tet_indices;

  bool is_valid() const {
    return num_degenerate_tets == 0 && num_inverted_tets == 0;
  }
};

/// Compute the signed volume of a tetrahedron given its 4 vertices.
/// Returns positive for correctly oriented tets, negative for inverted.
double tet_volume(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                  const Eigen::Vector3d &p2, const Eigen::Vector3d &p3);

/// Compute the volume of a tetrahedron from mesh element.
double tet_volume(const Mesh &mesh, const Element &tet);

/// Compute the aspect ratio of a tetrahedron (max edge / min edge).
double tet_aspect_ratio(const Mesh &mesh, const Element &tet);

/// Validate mesh quality, checking for degenerate and inverted elements.
MeshQualityReport validate_mesh(const Mesh &mesh);

} // namespace edgefem
