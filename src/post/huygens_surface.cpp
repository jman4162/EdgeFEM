#include "edgefem/post/huygens_surface.hpp"
#include "edgefem/edge_basis.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace edgefem {

namespace {

/// Make a sorted triple key for triangle-to-tet lookup.
/// Uses 21-bit packing: key = (a << 42) | (b << 21) | c where a < b < c.
uint64_t make_tri_key(int64_t a, int64_t b, int64_t c) {
  // Sort three values
  if (a > b)
    std::swap(a, b);
  if (b > c)
    std::swap(b, c);
  if (a > b)
    std::swap(a, b);
  return (static_cast<uint64_t>(a) << 42) | (static_cast<uint64_t>(b) << 21) |
         static_cast<uint64_t>(c);
}

} // namespace

HuygensSurfaceData extract_huygens_surface(const Mesh &mesh,
                                           const VecC &solution,
                                           int surface_tag, double omega,
                                           std::complex<double> mu_r) {
  const double mu0 = 4.0e-7 * M_PI;
  const std::complex<double> j_omega_mu =
      std::complex<double>(0, 1) * omega * mu0 * mu_r;

  // Step 1: Build lookup from sorted node triple -> tet index
  // Each tet has 4 faces; store face -> (tet_index, local_face_index)
  std::unordered_map<uint64_t, size_t> tri_to_tet;

  // Tet face connectivity: faces are opposite to vertices 0,1,2,3
  // Face 0 (opposite vertex 0): nodes 1,2,3
  // Face 1 (opposite vertex 1): nodes 0,2,3
  // Face 2 (opposite vertex 2): nodes 0,1,3
  // Face 3 (opposite vertex 3): nodes 0,1,2
  constexpr int face_nodes[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

  for (size_t t = 0; t < mesh.tets.size(); ++t) {
    const auto &tet = mesh.tets[t];
    for (int f = 0; f < 4; ++f) {
      uint64_t key =
          make_tri_key(tet.conn[face_nodes[f][0]], tet.conn[face_nodes[f][1]],
                       tet.conn[face_nodes[f][2]]);
      tri_to_tet[key] = t;
    }
  }

  HuygensSurfaceData data;

  for (const auto &tri : mesh.tris) {
    if (tri.phys != surface_tag)
      continue;

    // Triangle nodes
    int idx0 = mesh.nodeIndex.at(tri.conn[0]);
    int idx1 = mesh.nodeIndex.at(tri.conn[1]);
    int idx2 = mesh.nodeIndex.at(tri.conn[2]);

    Eigen::Vector3d v0 = mesh.nodes[idx0].xyz;
    Eigen::Vector3d v1 = mesh.nodes[idx1].xyz;
    Eigen::Vector3d v2 = mesh.nodes[idx2].xyz;

    // Centroid
    Eigen::Vector3d centroid = (v0 + v1 + v2) / 3.0;

    // Area and outward normal
    Eigen::Vector3d e1 = v1 - v0;
    Eigen::Vector3d e2 = v2 - v0;
    Eigen::Vector3d cross_vec = e1.cross(e2);
    double area = 0.5 * cross_vec.norm();
    Eigen::Vector3d normal = cross_vec.normalized();

    if (area < 1e-30)
      continue;

    // Find parent tet
    uint64_t key = make_tri_key(tri.conn[0], tri.conn[1], tri.conn[2]);
    auto it = tri_to_tet.find(key);
    if (it == tri_to_tet.end()) {
      continue; // Orphan triangle, skip
    }
    size_t tet_idx = it->second;
    const auto &tet = mesh.tets[tet_idx];

    // Get tet vertices
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int ni = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[ni].xyz;
    }

    // Ensure normal points outward (away from tet centroid)
    Eigen::Vector3d tet_centroid = (X[0] + X[1] + X[2] + X[3]) / 4.0;
    if (normal.dot(centroid - tet_centroid) < 0) {
      normal = -normal;
    }

    // Extract edge DOFs for this tet
    std::array<std::complex<double>, 6> edge_dofs;
    for (int e = 0; e < 6; ++e) {
      edge_dofs[e] = solution(tet.edges[e]);
    }

    // Evaluate E at centroid using Whitney interpolation
    Eigen::Vector3cd E =
        evaluate_edge_field(X, tet.edge_orient, edge_dofs, centroid);

    // Compute H = (1 / j*omega*mu) * curl(E)
    // For Whitney elements, curl is piecewise constant per tet
    Eigen::Matrix<double, 6, 3> curls = whitney_edge_curls(X);
    Eigen::Vector3cd curl_E = Eigen::Vector3cd::Zero();
    for (int e = 0; e < 6; ++e) {
      std::complex<double> coeff =
          edge_dofs[e] * static_cast<double>(tet.edge_orient[e]);
      curl_E += coeff * curls.row(e).transpose().cast<std::complex<double>>();
    }
    Eigen::Vector3cd H = curl_E / j_omega_mu;

    // Project to tangential components: X_tan = X - (X . n) * n
    Eigen::Vector3cd E_tan = E - (E.dot(normal.cast<std::complex<double>>())) *
                                     normal.cast<std::complex<double>>();
    Eigen::Vector3cd H_tan = H - (H.dot(normal.cast<std::complex<double>>())) *
                                     normal.cast<std::complex<double>>();

    data.r.push_back(centroid);
    data.n.push_back(normal);
    data.E_tan.push_back(E_tan);
    data.H_tan.push_back(H_tan);
    data.area.push_back(area);
  }

  if (data.r.empty()) {
    throw std::runtime_error(
        "extract_huygens_surface: no triangles found with surface_tag = " +
        std::to_string(surface_tag));
  }

  return data;
}

} // namespace edgefem
