#include "edgefem/ports/lumped_port.hpp"

#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace edgefem {

WavePort build_lumped_port(const Mesh &mesh, const LumpedPortConfig &config) {
  // Collect unique edges on the port surface
  std::unordered_set<int> edge_set;
  for (const auto &tri : mesh.tris) {
    if (tri.phys != config.surface_tag)
      continue;
    // Tri3 uses first 3 edges
    for (int e = 0; e < 3; ++e) {
      edge_set.insert(tri.edges[e]);
    }
  }

  if (edge_set.empty()) {
    throw std::runtime_error(
        "build_lumped_port: no triangles found with surface_tag = " +
        std::to_string(config.surface_tag));
  }

  // Build edge list and compute weights
  std::vector<int> edges(edge_set.begin(), edge_set.end());
  Eigen::VectorXcd weights(edges.size());

  Eigen::Vector3d e_dir = config.e_direction.normalized();

  for (size_t i = 0; i < edges.size(); ++i) {
    int ei = edges[i];
    const Edge &edge = mesh.edges[ei];

    // Get node positions
    int idx0 = mesh.nodeIndex.at(edge.n0);
    int idx1 = mesh.nodeIndex.at(edge.n1);
    Eigen::Vector3d p0 = mesh.nodes[idx0].xyz;
    Eigen::Vector3d p1 = mesh.nodes[idx1].xyz;

    // Edge vector (canonical direction: n0 < n1)
    Eigen::Vector3d edge_vec = p1 - p0;

    // Weight = projection of edge vector onto E-field direction
    double w = edge_vec.dot(e_dir);
    weights(i) = std::complex<double>(w, 0.0);
  }

  // Normalize weights for proper power coupling:
  // ||w||^2 should give sqrt(Z0) for consistent S-parameter extraction
  double w_norm = weights.norm();
  if (w_norm > 1e-15) {
    weights *= std::sqrt(config.z0) / w_norm;
  }

  WavePort port;
  port.surface_tag = config.surface_tag;
  port.edges = edges;
  port.weights = weights;
  port.mode.Z0 = config.z0;
  port.mode.beta = 0.0; // No propagation constant for lumped port
  port.mode.kc = 0.0;
  port.mode.fc = 0.0;

  return port;
}

} // namespace edgefem
