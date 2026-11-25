// calculate_sparams_eigenmode edge cases on WR-90 (examples/rect_waveguide.msh,
// 50 mm, Port1 at z = 0, Port2 at z = 50 mm).
//
// 1. Below the TE10 cutoff (6.557 GHz) the S-matrix is all NaN, not
//    uninitialized memory.
// 2. The port mode's sign (or phase) from the eigensolver is arbitrary; S21
//    must not depend on it. Negating one port's weights must leave S
//    unchanged.

#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iostream>

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

namespace {
struct Guide {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);
  RectWaveguidePort dims{0.02286, 0.01016};
};

std::vector<WavePort> ports_at(const Guide &g, double freq) {
  PortMode mode = solve_te10_mode(g.dims, freq);
  const double kc_sq = (M_PI / g.dims.a) * (M_PI / g.dims.a);
  return {build_wave_port_2d(g.mesh, 2, mode, g.bc.dirichlet_edges, kc_sq),
          build_wave_port_2d(g.mesh, 3, mode, g.bc.dirichlet_edges, kc_sq)};
}

Eigen::MatrixXcd sparams(const Guide &g, const std::vector<WavePort> &ports,
                         double freq) {
  MaxwellParams p;
  p.omega = 2.0 * M_PI * freq;
  return calculate_sparams_eigenmode(g.mesh, p, g.bc, ports);
}
} // namespace

int main() {
  Guide g;

  // 1. Below cutoff
  {
    const double f = 5e9;
    Eigen::MatrixXcd S = sparams(g, ports_at(g, 10e9), f);
    CHECK(S.rows() == 2 && S.cols() == 2);
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        CHECK(std::isnan(S(i, j).real()) && std::isnan(S(i, j).imag()));
  }

  // 2. Port-mode sign does not change S
  {
    const double f = 10e9;
    std::vector<WavePort> ports = ports_at(g, f);
    const Eigen::MatrixXcd S = sparams(g, ports, f);
    std::vector<WavePort> flipped2 = ports;
    flipped2[1].weights = -flipped2[1].weights;
    const Eigen::MatrixXcd S2 = sparams(g, flipped2, f);
    std::vector<WavePort> flipped1 = ports;
    flipped1[0].weights *= std::complex<double>(0.0, 1.0); // any phase
    const Eigen::MatrixXcd S1 = sparams(g, flipped1, f);
    std::cout << "S21 = " << S(1, 0) << ", port 2 negated: " << S2(1, 0)
              << ", port 1 times j: " << S1(1, 0) << "\n";
    CHECK((S - S2).cwiseAbs().maxCoeff() < 1e-9);
    CHECK((S - S1).cwiseAbs().maxCoeff() < 1e-9);
    // And the result is the physical one: phase(S21) = -beta L.
    const double k0 = 2.0 * M_PI * f / 299792458.0;
    const double kc = M_PI / g.dims.a;
    const double beta = std::sqrt(k0 * k0 - kc * kc);
    const double err = std::abs(
        std::arg(S(1, 0) * std::exp(std::complex<double>(0, beta * 0.05))));
    std::cout << "phase(S21) error vs -beta L: " << err * 180 / M_PI
              << " deg\n";
    CHECK(err < 10.0 * M_PI / 180.0);
  }

  // 3. One assembly and one factorization for all ports (direct solver),
  //    and the fields variant returns the same S with one field per port.
  {
    const double f = 10e9;
    std::vector<WavePort> ports = ports_at(g, f);
    MaxwellParams p;
    p.omega = 2.0 * M_PI * f;
    SolveOptions direct;
    direct.use_direct = true;
    const auto before = direct_factorization_count();
    SParamSolution sol =
        calculate_sparams_eigenmode_with_fields(g.mesh, p, g.bc, ports, direct);
    CHECK(direct_factorization_count() - before == 1);
    const Eigen::MatrixXcd S_plain =
        calculate_sparams_eigenmode(g.mesh, p, g.bc, ports);
    CHECK((sol.S - S_plain).cwiseAbs().maxCoeff() < 1e-8);
    CHECK(sol.x.size() == 2);
    for (const auto &x : sol.x) {
      CHECK(x.size() == static_cast<int>(g.mesh.edges.size()));
      CHECK(x.size() > 0 && x.norm() > 0.0);
      for (int e : g.bc.dirichlet_edges) // PEC edges carry no field
        if (x.size() > e)
          CHECK(x(e) == std::complex<double>(0.0));
    }
  }

  return edgefem_test::exit_code();
}
