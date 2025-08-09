#pragma once

#include <unordered_set>

#include "vectorem/mesh.hpp"

namespace vectorem {

struct BC {
  std::unordered_set<int> dirichlet_nodes;
  std::unordered_set<int> dirichlet_edges;
};

BC build_scalar_pec(const Mesh &mesh, int pec_tag);
BC build_edge_pec(const Mesh &mesh, int pec_tag);

} // namespace vectorem
