// Port-kind consistency of the S-parameter functions on WR-90.
//
// 1. Lumped and Modal ports use different S definitions (fill_sparam_column
//    picks the formula from the active port only), so a mix of kinds gives
//    a matrix that is neither one nor the other: the S-parameter functions
//    must refuse it.
// 2. The legacy diagonal port ABC (use_port_abc) must enter both the fast
//    (K/M) and slow frequency_sweep paths, so frequency_sweep at f equals
//    calculate_sparams at f.

#include "test_check.hpp"

#include <complex>
#include <iostream>
#include <stdexcept>

#include "edgefem/maxwell.hpp"
#include "edgefem/sweep.hpp"

using namespace edgefem;

namespace {
template <typename F> bool throws_invalid(F &&f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return true;
  }
  return false;
}
} // namespace

int main() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);
  const double f = 10e9;
  RectWaveguidePort dims{0.02286, 0.01016};
  PortMode mode = solve_te10_mode(dims, f);
  const double kc_sq = (M_PI / dims.a) * (M_PI / dims.a);
  WavePort modal1 =
      build_wave_port_2d(mesh, 2, mode, bc.dirichlet_edges, kc_sq);
  WavePort modal2 =
      build_wave_port_2d(mesh, 3, mode, bc.dirichlet_edges, kc_sq);
  WavePort lumped2 = modal2;
  lumped2.kind = PortKind::Lumped;

  MaxwellParams p;
  p.omega = 2.0 * M_PI * f;
  SolveOptions opts;
  opts.use_direct = true;

  // 1. Mixed kinds are refused.
  const std::vector<WavePort> mixed{modal1, lumped2};
  CHECK(throws_invalid([&] { calculate_sparams(mesh, p, bc, mixed, opts); }));
  CHECK(
      throws_invalid([&] { frequency_sweep(mesh, p, bc, mixed, {f}, opts); }));
  CHECK(throws_invalid(
      [&] { calculate_sparams_periodic(mesh, p, bc, PeriodicBC{}, mixed); }));

  // 2. use_port_abc enters both sweep paths.
  MaxwellParams p_abc = p;
  p_abc.use_port_abc = true;
  p_abc.port_abc_type = PortABCType::Beta;
  const std::vector<WavePort> modal{modal1, modal2};
  const Eigen::MatrixXcd S_direct =
      calculate_sparams(mesh, p_abc, bc, modal, opts);
  const SweepResult sweep = frequency_sweep(mesh, p_abc, bc, modal, {f}, opts);
  CHECK(sweep.S_matrices.size() == 1);
  if (sweep.S_matrices.size() == 1) {
    const double diff = (sweep.S_matrices[0] - S_direct).cwiseAbs().maxCoeff();
    std::cout << "use_port_abc: |S_sweep - S_calculate| = " << diff << "\n";
    CHECK(diff < 1e-9);
  }

  return edgefem_test::exit_code();
}
