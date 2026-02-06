#include "edgefem/mesh_quality.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>

namespace edgefem {

double tet_volume(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                  const Eigen::Vector3d &p2, const Eigen::Vector3d &p3) {
  // Volume = (1/6) * |det([p1-p0, p2-p0, p3-p0])|
  // Signed volume is positive for correctly oriented tets.
  Eigen::Vector3d v1 = p1 - p0;
  Eigen::Vector3d v2 = p2 - p0;
  Eigen::Vector3d v3 = p3 - p0;
  return v1.dot(v2.cross(v3)) / 6.0;
}

double tet_volume(const Mesh &mesh, const Element &tet) {
  const auto &n0 = mesh.nodes[mesh.nodeIndex.at(tet.conn[0])].xyz;
  const auto &n1 = mesh.nodes[mesh.nodeIndex.at(tet.conn[1])].xyz;
  const auto &n2 = mesh.nodes[mesh.nodeIndex.at(tet.conn[2])].xyz;
  const auto &n3 = mesh.nodes[mesh.nodeIndex.at(tet.conn[3])].xyz;
  return tet_volume(n0, n1, n2, n3);
}

double tet_aspect_ratio(const Mesh &mesh, const Element &tet) {
  const auto &n0 = mesh.nodes[mesh.nodeIndex.at(tet.conn[0])].xyz;
  const auto &n1 = mesh.nodes[mesh.nodeIndex.at(tet.conn[1])].xyz;
  const auto &n2 = mesh.nodes[mesh.nodeIndex.at(tet.conn[2])].xyz;
  const auto &n3 = mesh.nodes[mesh.nodeIndex.at(tet.conn[3])].xyz;

  // Compute all 6 edge lengths
  double edges[6];
  edges[0] = (n1 - n0).norm();
  edges[1] = (n2 - n0).norm();
  edges[2] = (n3 - n0).norm();
  edges[3] = (n2 - n1).norm();
  edges[4] = (n3 - n1).norm();
  edges[5] = (n3 - n2).norm();

  double min_edge = *std::min_element(edges, edges + 6);
  double max_edge = *std::max_element(edges, edges + 6);

  if (min_edge < std::numeric_limits<double>::epsilon()) {
    return std::numeric_limits<double>::infinity();
  }
  return max_edge / min_edge;
}

MeshQualityReport validate_mesh(const Mesh &mesh) {
  MeshQualityReport report;
  report.min_volume = std::numeric_limits<double>::max();
  report.max_aspect_ratio = 0.0;

  constexpr double volume_tol = 1e-15;

  for (size_t i = 0; i < mesh.tets.size(); ++i) {
    const auto &tet = mesh.tets[i];
    double vol = tet_volume(mesh, tet);
    double ar = tet_aspect_ratio(mesh, tet);

    if (vol < report.min_volume) {
      report.min_volume = vol;
    }
    if (ar > report.max_aspect_ratio) {
      report.max_aspect_ratio = ar;
    }

    if (std::abs(vol) < volume_tol) {
      report.num_degenerate_tets++;
      report.bad_tet_indices.push_back(static_cast<int>(i));
    } else if (vol < 0) {
      report.num_inverted_tets++;
      report.bad_tet_indices.push_back(static_cast<int>(i));
    }
  }

  return report;
}

} // namespace edgefem
