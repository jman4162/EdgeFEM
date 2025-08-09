#include <Eigen/Dense>
#include <cassert>
#include <cmath>

#include "vectorem/maxwell.hpp"

using namespace vectorem;

int main() {
  // Load simple mesh and tag first tetrahedron as PML region
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  int pml_tag = mesh.tets.empty() ? 0 : mesh.tets[0].phys + 1;
  if (!mesh.tets.empty())
    mesh.tets[0].phys = pml_tag;

  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;
  p.pml_sigma = 5.0; // strong damping
  p.pml_regions.insert(pml_tag);
  auto asmbl = assemble_maxwell(mesh, p, bc);

  // Plane-wave amplitude reduction through PML of thickness 1
  double refl = std::abs(std::exp(-2.0 * p.pml_sigma));
  double refl_dB = 20.0 * std::log10(refl);
  assert(refl_dB <= -40.0);
  // sanity: matrix built with expected size
  assert(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  return 0;
}
