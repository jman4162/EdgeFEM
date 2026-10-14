#include "test_check.hpp"
#include <array>

#include "edgefem/mesh.hpp"

using namespace edgefem;

int main() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  CHECK(!mesh.edges.empty());
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    const auto &e = mesh.edges[i];
    CHECK(e.n0 < e.n1);
    auto key = make_edge_key(e.n0, e.n1);
    CHECK(mesh.edgeIndex.at(key) == static_cast<int>(i));
  }
  const std::array<std::array<int, 2>, 6> tet_pairs = {
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2},
      std::array<int, 2>{0, 3}, std::array<int, 2>{1, 2},
      std::array<int, 2>{1, 3}, std::array<int, 2>{2, 3}};
  for (const auto &tet : mesh.tets) {
    for (int k = 0; k < 6; ++k) {
      auto a = tet.conn[tet_pairs[k][0]];
      auto b = tet.conn[tet_pairs[k][1]];
      int eid = tet.edges[k];
      int sign = tet.edge_orient[k];
      const auto &edge = mesh.edges[eid];
      CHECK(edge.n0 == std::min(a, b) && edge.n1 == std::max(a, b));
      int expected = (a < b) ? 1 : -1;
      CHECK(sign == expected);
    }
  }
  return edgefem_test::exit_code();
}
