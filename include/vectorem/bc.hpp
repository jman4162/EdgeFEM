#pragma once

#include <unordered_set>

#include "vectorem/mesh.hpp"

namespace vectorem {

struct BC {
  std::unordered_set<int> dirichlet_nodes;
};

BC build_scalar_pec(const Mesh &mesh, int pec_tag);

} // namespace vectorem
