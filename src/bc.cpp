#include "vectorem/bc.hpp"

#include <stdexcept>

namespace vectorem {

BC build_scalar_pec(const Mesh &mesh, int pec_tag) {
  BC bc;
  for (const auto &tri : mesh.tris) {
    if (tri.phys == pec_tag) {
      for (int k = 0; k < 3; ++k) {
        auto it = mesh.nodeIndex.find(tri.conn[k]);
        if (it != mesh.nodeIndex.end()) {
          bc.dirichlet_nodes.insert(it->second);
        }
      }
    }
  }
  return bc;
}

BC build_edge_pec(const Mesh &mesh, int pec_tag) {
  BC bc;
  for (const auto &tri : mesh.tris) {
    if (tri.phys == pec_tag) {
      for (int k = 0; k < 3; ++k) {
        bc.dirichlet_edges.insert(tri.edges[k]);
      }
    }
  }
  return bc;
}

} // namespace vectorem
