#include "edgefem/ports/lumped_port.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace edgefem {

// Compute w_i = ∫(N_i · e_dir) dS over port surface triangles.
// For a Whitney edge basis N_e on a triangle with nodes (v0, v1, v2),
// the edge (ni, nj) contributes:
//   ∫ (λ_i ∇λ_j - λ_j ∇λ_i) · e_dir dS
// On a flat triangle, ∇λ_i are constant and ∫ λ_i λ_j dS has closed forms:
//   ∫ λ_i² dS = A/6, ∫ λ_i λ_j dS = A/12 (i≠j)
// So ∫ N_e · e_dir dS = (1/6 - 1/12)(∇λ_j - ∇λ_i)·e_dir * sign
//                      = (A/12)(∇λ_j · e_dir - ∇λ_i · e_dir) * sign
// where ∇λ_k are the 2D surface gradients on the triangle.
static Eigen::VectorXcd compute_surface_integral_weights(
    const Mesh &mesh, int surface_tag, const Eigen::Vector3d &e_dir,
    const std::vector<int> &edges,
    const std::unordered_map<int, size_t> &edge_to_idx) {

  Eigen::VectorXcd weights = Eigen::VectorXcd::Zero(edges.size());

  for (const auto &tri : mesh.tris) {
    if (tri.phys != surface_tag)
      continue;

    // Get triangle vertex positions (Tri3 uses first 3 of conn)
    Eigen::Vector3d v[3];
    int node_ids[3];
    for (int k = 0; k < 3; ++k) {
      node_ids[k] = static_cast<int>(tri.conn[k]);
      int idx = mesh.nodeIndex.at(tri.conn[k]);
      v[k] = mesh.nodes[idx].xyz;
    }

    // Triangle area and normal
    Eigen::Vector3d e01 = v[1] - v[0];
    Eigen::Vector3d e02 = v[2] - v[0];
    Eigen::Vector3d normal = e01.cross(e02);
    double area2 = normal.norm(); // 2 * area
    if (area2 < 1e-30)
      continue;
    Eigen::Vector3d n_hat = normal.normalized();
    double area = area2 / 2.0;

    // 2D surface gradients of barycentric coordinates on the triangle.
    // For triangle (v0, v1, v2):
    //   ∇λ_0 = (n × (v2 - v1)) / (2A)
    //   ∇λ_1 = (n × (v0 - v2)) / (2A)
    //   ∇λ_2 = (n × (v1 - v0)) / (2A)
    Eigen::Vector3d edge_20 = v[2] - v[1];
    Eigen::Vector3d edge_01 = v[0] - v[2];
    Eigen::Vector3d edge_12 = v[1] - v[0];
    Eigen::Vector3d grad_lam[3];
    grad_lam[0] = n_hat.cross(edge_20) / area2;
    grad_lam[1] = n_hat.cross(edge_01) / area2;
    grad_lam[2] = n_hat.cross(edge_12) / area2;

    // Local edge table for Tri3: (n0_local, n1_local) for edges 0,1,2
    // Edge 0: (0,1), Edge 1: (0,2), Edge 2: (1,2)
    static const int tri_edge_nodes[3][2] = {{0, 1}, {0, 2}, {1, 2}};

    for (int le = 0; le < 3; ++le) {
      int ge = tri.edges[le];
      auto it = edge_to_idx.find(ge);
      if (it == edge_to_idx.end())
        continue;

      int li = tri_edge_nodes[le][0]; // local node index of edge start
      int lj = tri_edge_nodes[le][1]; // local node index of edge end

      // Whitney basis integral: ∫ N_e · e_dir dS = (A/12)(∇λ_j - ∇λ_i) · e_dir
      // But we also need orientation: the global edge is stored (n0 < n1).
      // Check if local orientation matches global.
      int orient = tri.edge_orient[le];

      // ∫(λ_i ∇λ_j - λ_j ∇λ_i) · e_dir dS
      // = (∇λ_j · e_dir) ∫ λ_i dS/1 - ... but with cross terms:
      // ∫ λ_i (∇λ_j · e_dir) dS = (∇λ_j · e_dir) * ∫ λ_i dS = (∇λ_j · e_dir) *
      // A/3 ∫ λ_j (∇λ_i · e_dir) dS = (∇λ_i · e_dir) * A/3 So integral =
      // (A/3)(∇λ_j · e_dir - ∇λ_i · e_dir)
      double integral =
          (area / 3.0) * ((grad_lam[lj] - grad_lam[li]).dot(e_dir));

      weights(it->second) += std::complex<double>(orient * integral, 0.0);
    }
  }

  return weights;
}

WavePort build_lumped_port(const Mesh &mesh, const LumpedPortConfig &config) {
  // Collect unique edges on the port surface
  std::unordered_set<int> edge_set;
  for (const auto &tri : mesh.tris) {
    if (tri.phys != config.surface_tag)
      continue;
    for (int e = 0; e < 3; ++e) {
      edge_set.insert(tri.edges[e]);
    }
  }

  if (edge_set.empty()) {
    throw std::runtime_error(
        "build_lumped_port: no triangles found with surface_tag = " +
        std::to_string(config.surface_tag));
  }

  // Build edge list and index map
  std::vector<int> edges(edge_set.begin(), edge_set.end());
  std::unordered_map<int, size_t> edge_to_idx;
  for (size_t i = 0; i < edges.size(); ++i) {
    edge_to_idx[edges[i]] = i;
  }

  Eigen::Vector3d e_dir = config.e_direction.normalized();
  Eigen::VectorXcd weights;

  if (config.weight_mode == LumpedPortWeightMode::SurfaceIntegral) {
    weights = compute_surface_integral_weights(mesh, config.surface_tag, e_dir,
                                               edges, edge_to_idx);
  } else {
    // Projection mode (original method)
    weights.resize(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
      int ei = edges[i];
      const Edge &edge = mesh.edges[ei];
      int idx0 = mesh.nodeIndex.at(edge.n0);
      int idx1 = mesh.nodeIndex.at(edge.n1);
      Eigen::Vector3d p0 = mesh.nodes[idx0].xyz;
      Eigen::Vector3d p1 = mesh.nodes[idx1].xyz;
      Eigen::Vector3d edge_vec = p1 - p0;
      weights(i) = std::complex<double>(edge_vec.dot(e_dir), 0.0);
    }
  }

  // Initial normalization: scale ||w|| to sqrt(Z0).
  // This is a rough starting point — normalize_port_weights() should be called
  // after assembly for accurate scaling via w^H A^{-1} w = Z0.
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
