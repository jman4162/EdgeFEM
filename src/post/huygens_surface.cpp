#include "edgefem/post/huygens_surface.hpp"
#include "edgefem/edge_basis.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

namespace edgefem {

namespace {

using FaceKey = std::array<std::int64_t, 3>;

FaceKey make_face_key(std::int64_t a, std::int64_t b, std::int64_t c) {
  FaceKey k{a, b, c};
  std::sort(k.begin(), k.end());
  return k;
}

struct FaceKeyHash {
  std::size_t operator()(const FaceKey &k) const {
    std::size_t h = std::hash<std::int64_t>()(k[0]);
    h ^= std::hash<std::int64_t>()(k[1]) + 0x9e3779b97f4a7c15ULL + (h << 6) +
         (h >> 2);
    h ^= std::hash<std::int64_t>()(k[2]) + 0x9e3779b97f4a7c15ULL + (h << 6) +
         (h >> 2);
    return h;
  }
};

// Up to two tets share a face.
struct Parents {
  std::size_t tet[2] = {0, 0};
  int count = 0;
};

// Projection onto the plane normal to n without conjugation (Eigen's dot()
// conjugates its first argument).
Eigen::Vector3cd tangential(const Eigen::Vector3cd &v,
                            const Eigen::Vector3d &n) {
  const std::complex<double> vn = v.x() * n.x() + v.y() * n.y() + v.z() * n.z();
  return v - vn * n.cast<std::complex<double>>();
}

} // namespace

HuygensSurfaceData
extract_huygens_surface(const Mesh &mesh, const VecC &solution, int surface_tag,
                        double omega, std::complex<double> mu_r,
                        std::optional<Eigen::Vector3d> interior_point) {
  const double mu0 = 4.0e-7 * M_PI;
  const std::complex<double> j_omega_mu =
      std::complex<double>(0, 1) * omega * mu0 * mu_r;

  // Face (sorted node triple) -> the one or two tets that share it.
  // Faces opposite vertices 0..3.
  constexpr int face_nodes[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  std::unordered_map<FaceKey, Parents, FaceKeyHash> face_to_tets;
  face_to_tets.reserve(2 * mesh.tets.size());
  for (std::size_t t = 0; t < mesh.tets.size(); ++t) {
    const auto &tet = mesh.tets[t];
    for (int f = 0; f < 4; ++f) {
      Parents &p = face_to_tets[make_face_key(tet.conn[face_nodes[f][0]],
                                              tet.conn[face_nodes[f][1]],
                                              tet.conn[face_nodes[f][2]])];
      if (p.count < 2)
        p.tet[p.count] = t;
      ++p.count;
    }
  }

  HuygensSurfaceData data;

  for (const auto &tri : mesh.tris) {
    if (tri.phys != surface_tag)
      continue;

    const Eigen::Vector3d v0 = mesh.nodes[mesh.nodeIndex.at(tri.conn[0])].xyz;
    const Eigen::Vector3d v1 = mesh.nodes[mesh.nodeIndex.at(tri.conn[1])].xyz;
    const Eigen::Vector3d v2 = mesh.nodes[mesh.nodeIndex.at(tri.conn[2])].xyz;
    const Eigen::Vector3d centroid = (v0 + v1 + v2) / 3.0;
    const Eigen::Vector3d cross_vec = (v1 - v0).cross(v2 - v0);
    const double area = 0.5 * cross_vec.norm();
    if (area < 1e-30)
      continue;
    Eigen::Vector3d normal = cross_vec.normalized();

    auto it =
        face_to_tets.find(make_face_key(tri.conn[0], tri.conn[1], tri.conn[2]));
    if (it == face_to_tets.end())
      continue; // orphan triangle
    const Parents &parents = it->second;
    const int n_par = std::min(parents.count, 2);

    // Orientation: away from interior_point when given; otherwise away from
    // the first parent tet, which is outward on a boundary face (one parent)
    // and arbitrary but deterministic on an interior face.
    if (interior_point) {
      if (normal.dot(centroid - *interior_point) < 0)
        normal = -normal;
    } else {
      const auto &tet = mesh.tets[parents.tet[0]];
      Eigen::Vector3d tc = Eigen::Vector3d::Zero();
      for (int i = 0; i < 4; ++i)
        tc += mesh.nodes[mesh.nodeIndex.at(tet.conn[i])].xyz / 4.0;
      if (normal.dot(centroid - tc) < 0)
        normal = -normal;
    }

    // E at the face centroid and H = -curl E / (j omega mu) (exp(+j omega t):
    // curl E = -j omega mu H), averaged over the parent tets. Tangential E is
    // continuous across the face; curl E is constant per tet, so the average
    // of the two sides is the better estimate of H on the face.
    Eigen::Vector3cd E = Eigen::Vector3cd::Zero();
    Eigen::Vector3cd H = Eigen::Vector3cd::Zero();
    for (int k = 0; k < n_par; ++k) {
      const auto &tet = mesh.tets[parents.tet[k]];
      std::array<Eigen::Vector3d, 4> X;
      for (int i = 0; i < 4; ++i)
        X[i] = mesh.nodes[mesh.nodeIndex.at(tet.conn[i])].xyz;
      std::array<std::complex<double>, 6> edge_dofs;
      for (int e = 0; e < 6; ++e)
        edge_dofs[e] = solution(tet.edges[e]);
      E += evaluate_edge_field(X, tet.edge_orient, edge_dofs, centroid);
      const Eigen::Matrix<double, 6, 3> curls = whitney_edge_curls(X);
      Eigen::Vector3cd curl_E = Eigen::Vector3cd::Zero();
      for (int e = 0; e < 6; ++e)
        curl_E += edge_dofs[e] * static_cast<double>(tet.edge_orient[e]) *
                  curls.row(e).transpose().cast<std::complex<double>>();
      H += -curl_E / j_omega_mu;
    }
    E /= static_cast<double>(n_par);
    H /= static_cast<double>(n_par);

    data.r.push_back(centroid);
    data.n.push_back(normal);
    data.E_tan.push_back(tangential(E, normal));
    data.H_tan.push_back(tangential(H, normal));
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
