#include "vectorem/ports/wave_port.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include "vectorem/edge_basis.hpp"

namespace vectorem {

namespace {
struct ArrayHash {
  std::size_t operator()(const std::array<std::int64_t, 2> &arr) const noexcept {
    return std::hash<std::int64_t>{}(arr[0]) ^ (std::hash<std::int64_t>{}(arr[1]) << 1);
  }
};

Eigen::Matrix<double, 3, 2>
compute_shape_gradients(const Mesh &mesh, const Element &tri) {
  const auto &n0 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[0]));
  const auto &n1 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[1]));
  const auto &n2 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[2]));

  Eigen::Matrix3d C;
  C << 1.0, n0.xyz.x(), n0.xyz.y(), 1.0, n1.xyz.x(), n1.xyz.y(), 1.0,
      n2.xyz.x(), n2.xyz.y();
  return C.inverse().block<3, 2>(0, 1);
}

Eigen::Vector2cd gradient_from_values(const Eigen::Matrix<double, 3, 2> &G,
                                      const Eigen::Vector3cd &values) {
  return G.transpose() * values;
}

Eigen::Matrix<double, 3, 2> to_shape_gradients(const Mesh &mesh,
                                               const Element &tri) {
  return compute_shape_gradients(mesh, tri);
}
} // namespace

PortSurfaceMesh extract_surface_mesh(const Mesh &volume_mesh, int surface_tag) {
  PortSurfaceMesh surface;
  std::unordered_map<std::int64_t, int> node_map;

  for (size_t tri_idx = 0; tri_idx < volume_mesh.tris.size(); ++tri_idx) {
    const auto &tri = volume_mesh.tris[tri_idx];
    if (tri.phys != surface_tag)
      continue;

    Element tri2;
    tri2.type = ElemType::Tri3;
    tri2.phys = tri.phys;

    for (int k = 0; k < 3; ++k) {
      auto node_id = tri.conn[k];
      auto it = node_map.find(node_id);
      if (it == node_map.end()) {
        const auto &node =
            volume_mesh.nodes.at(volume_mesh.nodeIndex.at(node_id));
        Node n;
        n.id = node.id;
        n.xyz = Eigen::Vector3d{node.xyz.x(), node.xyz.y(), 0.0};
        int new_idx = static_cast<int>(surface.mesh.nodes.size());
        surface.mesh.nodeIndex[n.id] = new_idx;
        surface.mesh.nodes.push_back(n);
        node_map[node_id] = new_idx;
      }
      tri2.conn[k] = tri.conn[k];
    }
    surface.mesh.tris.push_back(tri2);
    surface.volume_tri_indices.push_back(static_cast<int>(tri_idx));
  }

  // Build boundary lines from edges that occur once
  std::unordered_map<std::array<std::int64_t, 2>, int, ArrayHash> edge_counts;
  for (const auto &tri : surface.mesh.tris) {
    std::array<std::int64_t, 3> ids{tri.conn[0], tri.conn[1], tri.conn[2]};
    for (int e = 0; e < 3; ++e) {
      std::int64_t a = ids[e];
      std::int64_t b = ids[(e + 1) % 3];
      std::array<std::int64_t, 2> key{std::min(a, b), std::max(a, b)};
      edge_counts[key] += 1;
    }
  }

  for (const auto &kv : edge_counts) {
    if (kv.second == 1) {
      BoundaryLine bl;
      bl.n0 = kv.first[0];
      bl.n1 = kv.first[1];
      bl.phys = 1; // PEC boundary by default
      surface.mesh.boundary_lines.push_back(bl);
    }
  }

  return surface;
}

WavePort build_wave_port(const Mesh &volume_mesh, const PortSurfaceMesh &surface,
                         const PortMode &mode) {
  if (mode.field.size() !=
      static_cast<int>(surface.mesh.nodes.size())) {
    throw std::runtime_error(
        "Port mode field size does not match surface mesh nodes");
  }

  std::map<int, std::complex<double>> accum;
  std::map<int, int> counts;
  std::complex<double> normal_accum = 0.0;

  for (size_t tri_idx = 0; tri_idx < surface.mesh.tris.size(); ++tri_idx) {
    const auto &tri2d = surface.mesh.tris[tri_idx];
    const auto &tri3d =
        volume_mesh.tris.at(surface.volume_tri_indices.at(tri_idx));

    Eigen::Matrix<double, 3, 2> G = to_shape_gradients(surface.mesh, tri2d);
    Eigen::Vector3cd values;
    for (int k = 0; k < 3; ++k) {
      int idx = surface.mesh.nodeIndex.at(tri2d.conn[k]);
      values[k] = mode.field(idx);
    }
    Eigen::Vector2cd grad = gradient_from_values(G, values);

    std::complex<double> beta = mode.beta;
    double kc = mode.kc;
    if (kc == 0.0) {
      throw std::runtime_error("Port mode has zero cutoff wavenumber");
    }

    const std::complex<double> j(0.0, 1.0);

    Eigen::Matrix<std::complex<double>, 3, 1> Et;
    if (mode.pol == ModePolarization::TE) {
      // TE mode: Et = -jωμ/kc² × (ẑ × ∇t Hz)
      // ẑ × ∇t Hz = ẑ × (∂Hz/∂x x̂ + ∂Hz/∂y ŷ) = -∂Hz/∂y x̂ + ∂Hz/∂x ŷ
      // Wait, ẑ × x̂ = ŷ and ẑ × ŷ = -x̂
      // So: ẑ × (∂Hz/∂x x̂ + ∂Hz/∂y ŷ) = ∂Hz/∂x ŷ - ∂Hz/∂y x̂
      // Therefore: Et = -jωμ/kc² × (∂Hz/∂x ŷ - ∂Hz/∂y x̂)
      //              = jωμ/kc² × ∂Hz/∂y x̂ - jωμ/kc² × ∂Hz/∂x ŷ
      // So: Et_x = jωμ/kc² × ∂Hz/∂y, Et_y = -jωμ/kc² × ∂Hz/∂x
      auto factor = j * mode.omega * mode.mu / (kc * kc);
      Et << factor * grad(1), -factor * grad(0), 0.0;
    } else {
      auto factor = -beta / (kc * kc);
      Et << factor * grad(0), factor * grad(1), 0.0;
    }

    const auto &n0 =
        volume_mesh.nodes.at(volume_mesh.nodeIndex.at(tri3d.conn[0])).xyz;
    const auto &n1 =
        volume_mesh.nodes.at(volume_mesh.nodeIndex.at(tri3d.conn[1])).xyz;
    const auto &n2 =
        volume_mesh.nodes.at(volume_mesh.nodeIndex.at(tri3d.conn[2])).xyz;
    Eigen::Vector3d normal = (n1 - n0).cross(n2 - n0);
    normal_accum += normal.z();

    for (int e = 0; e < 3; ++e) {
      int edge_idx = tri3d.edges[e];
      const auto &edge = volume_mesh.edges.at(edge_idx);
      const auto &pa =
          volume_mesh.nodes.at(volume_mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &pb =
          volume_mesh.nodes.at(volume_mesh.nodeIndex.at(edge.n1)).xyz;
      Eigen::Vector3d edge_vec = pb - pa;
      std::complex<double> integral =
          Et(0) * edge_vec.x() + Et(1) * edge_vec.y() + Et(2) * edge_vec.z();
      accum[edge_idx] += integral;
      counts[edge_idx] += 1;
    }
  }

  double normal_sign = (std::real(normal_accum) >= 0.0) ? 1.0 : -1.0;

  WavePort port;
  port.surface_tag = surface.mesh.tris.empty() ? 0
                                                : surface.mesh.tris.front().phys;
  port.mode = mode;
  port.edges.reserve(accum.size());
  port.weights.resize(static_cast<int>(accum.size()));

  int idx = 0;
  for (const auto &kv : accum) {
    port.edges.push_back(kv.first);
  }
  std::sort(port.edges.begin(), port.edges.end());
  for (int edge_idx : port.edges) {
    std::complex<double> value = accum[edge_idx] / static_cast<double>(counts[edge_idx]);
    port.weights(idx++) = normal_sign * value;
  }

  return port;
}

void populate_te10_field(const PortSurfaceMesh &surface,
                         const RectWaveguidePort &port,
                         PortMode &mode) {
  const int num_nodes = static_cast<int>(surface.mesh.nodes.size());
  mode.field.resize(num_nodes);

  // Find the coordinate range of the surface
  double x_min = std::numeric_limits<double>::max();
  for (const auto &node : surface.mesh.nodes) {
    x_min = std::min(x_min, node.xyz.x());
  }

  // TE10 mode: Hz(x,y) = A * cos(π(x-x_min)/a)
  // The scalar potential from which Et is derived via Et = jωμ₀/kc² * ∇t Hz
  for (int i = 0; i < num_nodes; ++i) {
    const auto &node = surface.mesh.nodes[i];
    double x = node.xyz.x() - x_min;
    double Hz = std::cos(M_PI * x / port.a);
    mode.field(i) = Hz;
  }

  // Normalize for unit power using mode parameters.
  // For TE10 mode (Pozar convention):
  //   Hz = A cos(πx/a)
  //   Ey = (jωμ/kc²)(π/a) A sin(πx/a)
  //   Hx = (jβ/kc²)(π/a) A sin(πx/a)
  //
  // Time-averaged power flow:
  // P = ∫∫ ½Re{Ey Hx*} dxdy
  //   = ½ (ωμβπ²/a²kc⁴) A² ∫sin²(πx/a)dx ∫dy
  //   = ½ (ωμβπ²/a²kc⁴) A² (a/2) b
  //   = ωμβπ²bA² / (4a·kc⁴)
  //
  // Using kc = π/a, so kc⁴ = π⁴/a⁴:
  // P = ωμβπ²bA² · a⁴ / (4a·π⁴) = ωμβbA²a³ / (4π²)
  //
  // For unit power: A² = 4π² / (ωμβba³)
  // Equivalently:   A² = 4·kc⁴·a / (ωμβπ²b)
  //
  // Note: mode.mu is the material μ = μ₀*μ_r (we use real part for scaling)
  double omega = mode.omega;
  double mu = std::real(mode.mu);
  double beta = std::real(mode.beta);
  double kc = mode.kc;
  double pi_sq = M_PI * M_PI;

  double A_sq = 4.0 * std::pow(kc, 4) * port.a / (omega * mu * beta * pi_sq * port.b);
  double scale = std::sqrt(A_sq);

  mode.field *= scale;
}

WavePort build_wave_port_from_eigenvector(const Mesh &volume_mesh,
                                           const PortSurfaceMesh &surface,
                                           const Eigen::VectorXd &eigenvector,
                                           const PortMode &mode,
                                           const std::unordered_set<int> &pec_edges) {
  WavePort port;
  port.surface_tag = surface.mesh.tris.empty() ? 0 : surface.mesh.tris.front().phys;
  port.mode = mode;

  // Collect all edges from surface triangles
  std::unordered_set<int> port_edge_set;
  for (size_t tri_idx = 0; tri_idx < surface.mesh.tris.size(); ++tri_idx) {
    const auto &tri3d = volume_mesh.tris.at(surface.volume_tri_indices.at(tri_idx));
    for (int e = 0; e < 3; ++e) {
      port_edge_set.insert(tri3d.edges[e]);
    }
  }

  // Sort edges for consistent ordering
  port.edges.assign(port_edge_set.begin(), port_edge_set.end());
  std::sort(port.edges.begin(), port.edges.end());

  // Extract eigenvector values at port edges
  port.weights.resize(port.edges.size());
  double norm_sq = 0.0;
  int free_count = 0;

  for (size_t i = 0; i < port.edges.size(); ++i) {
    int e = port.edges[i];
    if (pec_edges.count(e)) {
      port.weights(i) = std::complex<double>(0.0, 0.0);
    } else {
      // Convert real eigenvector to complex (multiply by j to match TE mode phase)
      port.weights(i) = std::complex<double>(0.0, eigenvector(e));
      norm_sq += eigenvector(e) * eigenvector(e);
      free_count++;
    }
  }

  // Normalize so ||w||² = sqrt(Z0) for proper S-parameter extraction
  // This ensures V = w^H * x = sqrt(Z0) for unit-amplitude mode
  if (norm_sq > 1e-15 && std::real(mode.Z0) > 1e-15) {
    double target_norm_sq = std::sqrt(std::real(mode.Z0));
    double scale_factor = std::sqrt(target_norm_sq / norm_sq);
    port.weights *= scale_factor;
  }

  return port;
}

Eigen::VectorXd compute_te_eigenvector(const Mesh &mesh,
                                        const std::unordered_set<int> &pec_edges,
                                        double target_kc_sq) {
  const int n = static_cast<int>(mesh.edges.size());

  // Build mapping from original to free DOF indices
  std::vector<int> free_to_orig;
  std::vector<int> orig_to_free(n, -1);
  for (int i = 0; i < n; ++i) {
    if (!pec_edges.count(i)) {
      orig_to_free[i] = static_cast<int>(free_to_orig.size());
      free_to_orig.push_back(i);
    }
  }
  const int n_free = static_cast<int>(free_to_orig.size());

  if (n_free == 0) {
    return Eigen::VectorXd::Zero(n);
  }

  // Assemble K and M matrices for free DOFs
  // Use very low frequency to get K, and difference method to get M
  constexpr double c0 = 299792458.0;
  constexpr double mu0 = 4.0 * M_PI * 1e-7;

  // Assemble at low frequency (essentially K matrix)
  Eigen::MatrixXd K_free = Eigen::MatrixXd::Zero(n_free, n_free);
  Eigen::MatrixXd M_free = Eigen::MatrixXd::Zero(n_free, n_free);

  // Build element matrices and assemble
  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }

    auto Kloc = whitney_curl_curl_matrix(X);
    auto Mloc = whitney_mass_matrix(X);

    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int fi = orig_to_free[gi];
      if (fi < 0) continue;  // PEC edge
      int si = tet.edge_orient[i];

      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int fj = orig_to_free[gj];
        if (fj < 0) continue;  // PEC edge
        int sj = tet.edge_orient[j];

        double sign = static_cast<double>(si * sj);
        K_free(fi, fj) += sign * Kloc(i, j);
        M_free(fi, fj) += sign * Mloc(i, j);
      }
    }
  }

  // Solve generalized eigenvalue problem K*v = lambda*M*v
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);

  if (ges.info() != Eigen::Success) {
    // Return zero vector on failure
    return Eigen::VectorXd::Zero(n);
  }

  const auto &eigenvalues = ges.eigenvalues();
  const auto &eigenvectors = ges.eigenvectors();

  // Find eigenvalue closest to target_kc_sq (skip near-zero eigenvalues)
  int best_idx = -1;
  double min_diff = std::numeric_limits<double>::max();

  for (int i = 0; i < n_free; ++i) {
    if (eigenvalues(i) > 100.0) {  // Skip null space
      double diff = std::abs(eigenvalues(i) - target_kc_sq);
      if (diff < min_diff) {
        min_diff = diff;
        best_idx = i;
      }
    }
  }

  if (best_idx < 0) {
    return Eigen::VectorXd::Zero(n);
  }

  // Expand eigenvector to full size (zero for PEC edges)
  Eigen::VectorXd v_full = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n_free; ++i) {
    v_full(free_to_orig[i]) = eigenvectors(i, best_idx);
  }

  return v_full;
}

} // namespace vectorem

