#include <cassert>
#include <Eigen/Dense>

#include "vectorem/maxwell.hpp"

using namespace vectorem;

int main() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;
  auto asmbl = assemble_maxwell(mesh, p, bc);
  assert(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  Eigen::MatrixXcd A = Eigen::MatrixXcd(asmbl.A);
  assert(A.isApprox(A.transpose(), 1e-9));
  return 0;
}
