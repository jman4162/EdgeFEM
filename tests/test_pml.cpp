#include <Eigen/Dense>
#include <algorithm>
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
  PMLRegionSpec spec;
  spec.sigma_max = Eigen::Vector3d::Constant(5.0);
  spec.thickness = Eigen::Vector3d::Constant(1.0);
  spec.grading_order = 2.0;
  p.pml_tensor_regions.emplace(pml_tag, spec);
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1); // No ports

  assert(!asmbl.diagnostics.empty());
  const auto &diag = asmbl.diagnostics.front();
  for (int axis = 0; axis < 3; ++axis) {
    double refl = diag.reflection_est[axis];
    double refl_dB = 20.0 * std::log10(std::max(refl, 1e-12));
    // Expect at least ~40 dB attenuation per axis
    assert(refl_dB <= -40.0);
  }
  // sanity: matrix built with expected size
  assert(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  return 0;
}
