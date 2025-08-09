#include "vectorem/mesh.hpp"

#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace vectorem {
namespace {
void build_edges(Mesh &mesh);
}

Mesh load_gmsh_v2(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open mesh file: " + path);
  }

  Mesh mesh;
  std::string line;
  while (std::getline(in, line)) {
    if (line == "$Nodes") {
      int n;
      in >> n;
      for (int i = 0; i < n; ++i) {
        std::int64_t id;
        double x, y, z;
        in >> id >> x >> y >> z;
        mesh.nodeIndex[id] = static_cast<int>(mesh.nodes.size());
        mesh.nodes.push_back(Node{id, Eigen::Vector3d{x, y, z}});
      }
      std::getline(in, line); // consume endline
      std::getline(in, line); // "$EndNodes"
    } else if (line == "$Elements") {
      int m;
      in >> m;
      for (int i = 0; i < m; ++i) {
        std::int64_t id;
        int type, num_tags;
        in >> id >> type >> num_tags;
        Element e;
        e.id = id;
        e.type = static_cast<ElemType>(type);
        for (int t = 0; t < num_tags; ++t) {
          int tag;
          in >> tag;
          if (t == 0)
            e.phys = tag;
        }
        if (type == 2) { // Tri3
          for (int k = 0; k < 3; ++k)
            in >> e.conn[k];
          mesh.tris.push_back(e);
        } else if (type == 4) { // Tet4
          for (int k = 0; k < 4; ++k)
            in >> e.conn[k];
          mesh.tets.push_back(e);
        } else {
          std::string dummy;
          std::getline(in, dummy);
        }
      }
      std::getline(in, line); // consume endline
      std::getline(in, line); // "$EndElements"
    }
  }
  build_edges(mesh);
  return mesh;
}

namespace {
void build_edges(Mesh &mesh) {
  const std::array<std::array<int, 2>, 6> tet_pairs = {
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2},
      std::array<int, 2>{0, 3}, std::array<int, 2>{1, 2},
      std::array<int, 2>{1, 3}, std::array<int, 2>{2, 3}};
  const std::array<std::array<int, 2>, 3> tri_pairs = {
      std::array<int, 2>{0, 1}, std::array<int, 2>{1, 2},
      std::array<int, 2>{2, 0}};
  for (auto &tet : mesh.tets) {
    for (int e = 0; e < 6; ++e) {
      auto a = tet.conn[tet_pairs[e][0]];
      auto b = tet.conn[tet_pairs[e][1]];
      int sign = (a < b) ? 1 : -1;
      std::uint64_t key = make_edge_key(a, b);
      auto it = mesh.edgeIndex.find(key);
      if (it == mesh.edgeIndex.end()) {
        int idx = static_cast<int>(mesh.edges.size());
        mesh.edgeIndex[key] = idx;
        mesh.edges.push_back(Edge{std::min(a, b), std::max(a, b)});
        it = mesh.edgeIndex.find(key);
      }
      tet.edges[e] = it->second;
      tet.edge_orient[e] = sign;
    }
  }
  for (auto &tri : mesh.tris) {
    for (int e = 0; e < 3; ++e) {
      auto a = tri.conn[tri_pairs[e][0]];
      auto b = tri.conn[tri_pairs[e][1]];
      int sign = (a < b) ? 1 : -1;
      std::uint64_t key = make_edge_key(a, b);
      auto it = mesh.edgeIndex.find(key);
      if (it == mesh.edgeIndex.end()) {
        int idx = static_cast<int>(mesh.edges.size());
        mesh.edgeIndex[key] = idx;
        mesh.edges.push_back(Edge{std::min(a, b), std::max(a, b)});
        it = mesh.edgeIndex.find(key);
      }
      tri.edges[e] = it->second;
      tri.edge_orient[e] = sign;
    }
  }
}
} // namespace
} // namespace vectorem
