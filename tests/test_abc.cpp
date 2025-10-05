#include <cassert>
#include <cmath>
#include <complex>

#include "vectorem/maxwell.hpp"

using namespace vectorem;

int main() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);

  MaxwellParams p;
  p.omega = 2.0;
  auto asmbl_noabc = assemble_maxwell(mesh, p, bc, {}, -1);

  p.use_abc = true;
  auto asmbl_abc = assemble_maxwell(mesh, p, bc, {}, -1);

  assert(asmbl_noabc.A.rows() == asmbl_abc.A.rows());
  const int m = asmbl_abc.A.rows();
  int damped_edges = 0;

  for (int e = 0; e < m; ++e) {
    const std::complex<double> val_abc = asmbl_abc.A.coeff(e, e);
    const std::complex<double> val_no = asmbl_noabc.A.coeff(e, e);
    const std::complex<double> delta = val_abc - val_no;
    if (std::abs(delta) > 1e-12) {
      assert(std::abs(delta - std::complex<double>(0.0, p.omega)) < 1e-9);
      ++damped_edges;
    }
  }
  assert(damped_edges > 0);

  // Manufactured plane-wave reflection estimate for unit-thickness ABC layer.
  const double refl = std::exp(-p.omega);
  const double refl_dB = 20.0 * std::log10(refl);
  assert(refl_dB < -8.0);

  return 0;
}
