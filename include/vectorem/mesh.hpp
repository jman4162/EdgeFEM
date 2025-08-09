#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace vectorem {

struct Node {
  std::int64_t id;
  Eigen::Vector3d xyz;
};

struct Edge {
  std::int64_t n0;
  std::int64_t n1;
};

enum class ElemType : int {
  Tri3 = 2, // surface triangle
  Tet4 = 4  // volume tetrahedron
};

struct Element {
  std::int64_t id;
  ElemType type;
  std::array<std::int64_t, 4> conn{}; // Tet4 uses first 4 entries
  int phys = 0;                       // physical tag (BC/material)
  std::array<int, 6> edges{};         // global edge indices
  std::array<int, 6> edge_orient{};   // +1 if same orientation as global
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<Element> tets;
  std::vector<Element> tris; // boundary faces
  // maps for quick lookup
  std::vector<Edge> edges;                          // global edges (n0 < n1)
  std::unordered_map<std::uint64_t, int> edgeIndex; // key -> edge index
  std::unordered_map<std::int64_t, int> nodeIndex;
};

Mesh load_gmsh_v2(const std::string &path);

inline std::uint64_t make_edge_key(std::int64_t a, std::int64_t b) {
  if (a > b)
    std::swap(a, b);
  return (static_cast<std::uint64_t>(a) << 32) ^ static_cast<std::uint64_t>(b);
}

} // namespace vectorem
