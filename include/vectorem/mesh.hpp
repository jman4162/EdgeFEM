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

enum class ElemType : int {
  Tri3 = 2, // surface triangle
  Tet4 = 4  // volume tetrahedron
};

struct Element {
  std::int64_t id;
  ElemType type;
  std::array<std::int64_t, 4> conn{}; // Tet4 uses first 4 entries
  int phys = 0;                       // physical tag (BC/material)
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<Element> tets;
  std::vector<Element> tris; // boundary faces
  // maps for quick lookup
  std::unordered_map<std::int64_t, int> nodeIndex;
};

Mesh load_gmsh_v2(const std::string &path);

} // namespace vectorem
