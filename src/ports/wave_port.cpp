#include "vectorem/ports/wave_port.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <map>
#include <stdexcept>
#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/Geometry>

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
      auto factor = j * mode.omega * mode.mu / (kc * kc);
      Et << -factor * grad(1), factor * grad(0), 0.0;
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

} // namespace vectorem

