// PMLDiagnostic::reflection_est against the closed form for a graded layer.
// With s = 1 - j sigma / omega a normally incident plane wave decays as
// exp(-int sigma dz / c0), so a layer sigma = sigma_max (z / d)^n reflects
//   R = exp(-2 sigma_max d / ((n + 1) c0))
// after the round trip (independent of frequency). tests/test_pml_waveguide
// checks the absorption in an actual solve.

#include "test_check.hpp"
#include <cmath>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

int main() {
  constexpr double c0 = 299792458.0;
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  CHECK(!mesh.tets.empty());
  const int pml_tag = mesh.tets[0].phys + 1;
  mesh.tets[0].phys = pml_tag;

  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 2.0 * M_PI * 10e9;

  // sigma_max chosen for R = 0.01 (-40 dB) on x; y half of it; z off.
  const double d = 0.02;
  const double n = 2.0;
  const double sigma_40dB = std::log(100.0) * (n + 1.0) * c0 / (2.0 * d);
  PMLRegionSpec spec;
  spec.sigma_max = Eigen::Vector3d(sigma_40dB, 0.5 * sigma_40dB, 0.0);
  spec.thickness = Eigen::Vector3d::Constant(d);
  spec.grading_order = n;
  p.pml_tensor_regions.emplace(pml_tag, spec);
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);

  CHECK(asmbl.diagnostics.size() == 1);
  const auto &diag = asmbl.diagnostics.front();
  CHECK(diag.region_tag == pml_tag);
  CHECK_NEAR(diag.reflection_est[0], 0.01, 1e-12);
  CHECK_NEAR(diag.reflection_est[1], 0.1, 1e-12); // sqrt(0.01)
  CHECK(diag.reflection_est[2] == 1.0);           // no loss on z

  // The estimate does not depend on frequency.
  p.omega = 2.0 * M_PI * 1e9;
  auto asmbl_1ghz = assemble_maxwell(mesh, p, bc, {}, -1);
  CHECK_NEAR(asmbl_1ghz.diagnostics.front().reflection_est[0], 0.01, 1e-12);

  CHECK(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  return edgefem_test::exit_code();
}
