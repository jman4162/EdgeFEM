#include "edgefem/periodic.hpp"

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace edgefem {

namespace {

// Compute edge centroid from mesh
Eigen::Vector3d edge_centroid(const Mesh &mesh, int edge_idx) {
  const Edge &e = mesh.edges[edge_idx];
  const Eigen::Vector3d &p0 = mesh.nodes[mesh.nodeIndex.at(e.n0)].xyz;
  const Eigen::Vector3d &p1 = mesh.nodes[mesh.nodeIndex.at(e.n1)].xyz;
  return 0.5 * (p0 + p1);
}

// Get all edges on a surface with given physical tag
std::vector<int> get_surface_edges(const Mesh &mesh, int surface_tag) {
  std::unordered_set<int> edge_set;

  for (const auto &tri : mesh.tris) {
    if (tri.phys == surface_tag) {
      // Triangle has 3 edges (first 3 in the edges array)
      for (int i = 0; i < 3; ++i) {
        edge_set.insert(tri.edges[i]);
      }
    }
  }

  return std::vector<int>(edge_set.begin(), edge_set.end());
}

// Compute edge orientation on a triangle
int edge_orientation_on_tri(const Element &tri, int edge_idx) {
  for (int i = 0; i < 3; ++i) {
    if (tri.edges[i] == edge_idx) {
      return tri.edge_orient[i];
    }
  }
  return 0; // Edge not found on this triangle
}

} // namespace

PeriodicBC build_periodic_pairs(const Mesh &mesh, int master_tag, int slave_tag,
                                const Eigen::Vector3d &period_vector,
                                double tolerance) {
  PeriodicBC pbc;
  pbc.period_vector = period_vector;
  pbc.phase_shift = std::complex<double>(1.0, 0.0); // Default: no phase shift

  // Get edges on master and slave surfaces
  std::vector<int> master_edges = get_surface_edges(mesh, master_tag);
  std::vector<int> slave_edges = get_surface_edges(mesh, slave_tag);

  if (master_edges.empty()) {
    throw std::runtime_error("No edges found on master surface with tag " +
                             std::to_string(master_tag));
  }
  if (slave_edges.empty()) {
    throw std::runtime_error("No edges found on slave surface with tag " +
                             std::to_string(slave_tag));
  }

  // Build spatial index for slave edges
  // Map from quantized position to edge index
  double scale = 1.0 / tolerance;
  std::unordered_map<std::int64_t, int> slave_map;

  auto quantize = [scale](const Eigen::Vector3d &p) -> std::int64_t {
    // Simple spatial hash
    int64_t ix = static_cast<int64_t>(std::round(p.x() * scale));
    int64_t iy = static_cast<int64_t>(std::round(p.y() * scale));
    int64_t iz = static_cast<int64_t>(std::round(p.z() * scale));
    // Mix the coordinates
    return (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791);
  };

  for (int slave_idx : slave_edges) {
    Eigen::Vector3d centroid = edge_centroid(mesh, slave_idx);
    slave_map[quantize(centroid)] = slave_idx;
  }

  // Match master edges to slave edges
  std::unordered_set<int> matched_slaves;

  for (int master_idx : master_edges) {
    Eigen::Vector3d master_centroid = edge_centroid(mesh, master_idx);
    Eigen::Vector3d expected_slave = master_centroid + period_vector;

    // Look up in spatial hash
    int64_t hash = quantize(expected_slave);
    auto it = slave_map.find(hash);

    int matched_slave = -1;

    if (it != slave_map.end()) {
      // Verify the match is within tolerance
      Eigen::Vector3d slave_centroid = edge_centroid(mesh, it->second);
      if ((slave_centroid - expected_slave).norm() < tolerance) {
        matched_slave = it->second;
      }
    }

    // If hash lookup failed, do brute force search (handles hash collisions)
    if (matched_slave < 0) {
      for (int slave_idx : slave_edges) {
        Eigen::Vector3d slave_centroid = edge_centroid(mesh, slave_idx);
        if ((slave_centroid - expected_slave).norm() < tolerance) {
          matched_slave = slave_idx;
          break;
        }
      }
    }

    if (matched_slave < 0) {
      throw std::runtime_error(
          "Could not find matching slave edge for master edge " +
          std::to_string(master_idx));
    }

    if (matched_slaves.count(matched_slave)) {
      throw std::runtime_error("Slave edge " + std::to_string(matched_slave) +
                               " matched to multiple master edges");
    }
    matched_slaves.insert(matched_slave);

    // Determine orientations
    int master_orient = 0;
    int slave_orient = 0;

    for (const auto &tri : mesh.tris) {
      if (tri.phys == master_tag) {
        int o = edge_orientation_on_tri(tri, master_idx);
        if (o != 0)
          master_orient = o;
      }
      if (tri.phys == slave_tag) {
        int o = edge_orientation_on_tri(tri, matched_slave);
        if (o != 0)
          slave_orient = o;
      }
    }

    // Default to +1 if not found (edge may be on boundary of surface)
    if (master_orient == 0)
      master_orient = 1;
    if (slave_orient == 0)
      slave_orient = 1;

    PeriodicPair pair;
    pair.master_edge = master_idx;
    pair.slave_edge = matched_slave;
    pair.master_orient = master_orient;
    pair.slave_orient = slave_orient;
    pair.translation = period_vector;
    pbc.pairs.push_back(pair);
  }

  return pbc;
}

bool validate_periodic_bc(const Mesh &mesh, const PeriodicBC &pbc) {
  if (pbc.pairs.empty()) {
    return false;
  }

  // Check that all pairs have valid edge indices
  int num_edges = static_cast<int>(mesh.edges.size());
  for (const auto &pair : pbc.pairs) {
    if (pair.master_edge < 0 || pair.master_edge >= num_edges) {
      return false;
    }
    if (pair.slave_edge < 0 || pair.slave_edge >= num_edges) {
      return false;
    }
    if (pair.master_orient != 1 && pair.master_orient != -1) {
      return false;
    }
    if (pair.slave_orient != 1 && pair.slave_orient != -1) {
      return false;
    }
  }

  // Check that no edge appears in multiple pairs
  std::unordered_set<int> master_set, slave_set;
  for (const auto &pair : pbc.pairs) {
    if (master_set.count(pair.master_edge)) {
      return false;
    }
    if (slave_set.count(pair.slave_edge)) {
      return false;
    }
    master_set.insert(pair.master_edge);
    slave_set.insert(pair.slave_edge);
  }

  return true;
}

void set_floquet_phase(PeriodicBC &pbc, const Eigen::Vector2d &k_transverse) {
  // Phase = exp(j * (kx*Lx + ky*Ly))
  // Assumes period_vector is in x-y plane primarily
  double phase_arg = k_transverse.x() * pbc.period_vector.x() +
                     k_transverse.y() * pbc.period_vector.y();
  pbc.phase_shift = std::complex<double>(std::cos(phase_arg), std::sin(phase_arg));
}

std::complex<double> floquet_phase_from_angle(
    const Eigen::Vector3d &period_vector, double theta, double phi, double k0) {
  // Incident wave vector: k = k0 * [sin(theta)cos(phi), sin(theta)sin(phi),
  // -cos(theta)]
  Eigen::Vector3d k_inc;
  k_inc.x() = k0 * std::sin(theta) * std::cos(phi);
  k_inc.y() = k0 * std::sin(theta) * std::sin(phi);
  k_inc.z() = -k0 * std::cos(theta); // Negative z for downward propagation

  double phase_arg = k_inc.dot(period_vector);
  return std::complex<double>(std::cos(phase_arg), std::sin(phase_arg));
}

int count_surface_edges(const Mesh &mesh, int surface_tag) {
  std::vector<int> edges = get_surface_edges(mesh, surface_tag);
  return static_cast<int>(edges.size());
}

} // namespace edgefem
