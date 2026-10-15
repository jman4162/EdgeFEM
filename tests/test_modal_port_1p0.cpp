// Modal ports keep the 1.0.0 port terms (load w w^H / Z0, source
// 2 w / sqrt(Z0), S(j, i) = V_j / sqrt(Z0_i) - delta_ij); only ports from
// build_lumped_port get the gap-port terms. This test pins the modal
// S-parameters on examples/rect_waveguide.msh (WR-90, 50 mm) to the values
// EdgeFEM 1.0.0 returned on the same mesh, through calculate_sparams and
// both frequency_sweep paths.
//
// The stored values are not physical: a matched lossless guide has
// S11 = 0, |S21| = 1, and 1.0.0 gives |S11| = 2.41 at 8 GHz with analytic
// TE10 weights and S11 = -1.0000 with 2D port-face modes. Modal-port
// S-parameters are not validated in this release; the test only guards
// against an unannounced change. The validated modal path is
// calculate_sparams_eigenmode (test_eigenmode_sparams).
//
// Reference values: pyedgefem 1.0.0 (Python 3.12 wheel built in August
// 2026), use_direct = true, run on 2026-10-15. The branch reproduced them
// to 1e-13 (iterative default solver: 1e-12).

#include "edgefem/maxwell.hpp"
#include "edgefem/sweep.hpp"
#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using namespace edgefem;
using cd = std::complex<double>;

namespace {

const std::vector<double> kFreqs = {8e9, 10e9, 12e9};

// 1.0.0: {S11, S21} per frequency, analytic TE10 weights (build_wave_port).
const cd kTe10[3][2] = {
    {{-2.408954280895234, 0.0}, {0.006938625438499754, 0.0}},
    {{-1.62805583153405, 0.0}, {-0.0024527851536667866, 0.0}},
    {{-1.3296373741166032, 0.0}, {0.0013598394437911014, 0.0}},
};

// 1.0.0: {S11, S21} per frequency, 2D port-face modes (build_wave_port_2d).
const cd kPort2d[3][2] = {
    {{-0.9999832537796255, 0.0}, {0.00010575256402561964, 0.0}},
    {{-0.999983560978861, 0.0}, {-8.548181105334578e-05, 0.0}},
    {{-1.000018303705424, 0.0}, {6.808798103803714e-05, 0.0}},
};

constexpr double kTol = 1e-9;

std::vector<WavePort> te10_ports(const Mesh &mesh,
                                 const RectWaveguidePort &dims, double f) {
  std::vector<WavePort> ports;
  for (int tag : {2, 3}) {
    PortSurfaceMesh surf = extract_surface_mesh(mesh, tag);
    PortMode mode = solve_te10_mode(dims, f);
    populate_te10_field(surf, dims, mode);
    ports.push_back(build_wave_port(mesh, surf, mode));
  }
  return ports;
}

void check_column(const Eigen::MatrixXcd &S, const cd (&ref)[2],
                  const char *what, double f) {
  std::cout << what << " " << f / 1e9 << " GHz: S11=" << S(0, 0)
            << " S21=" << S(1, 0) << std::endl;
  CHECK_NEAR(S(0, 0), ref[0], kTol);
  CHECK_NEAR(S(1, 0), ref[1], kTol);
}

} // namespace

int main() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);
  RectWaveguidePort dims{0.02286, 0.01016};
  SolveOptions opts;
  opts.use_direct = true;

  // Every builder except build_lumped_port gives a Modal port.
  {
    const auto ports = te10_ports(mesh, dims, 10e9);
    CHECK(ports[0].kind == PortKind::Modal);
    CHECK(has_modal_ports(ports));
  }

  for (size_t i = 0; i < kFreqs.size(); ++i) {
    const double f = kFreqs[i];
    MaxwellParams p;
    p.omega = 2 * M_PI * f;
    check_column(
        calculate_sparams(mesh, p, bc, te10_ports(mesh, dims, f), opts),
        kTe10[i], "calculate_sparams TE10", f);

    const PortMode mode = solve_te10_mode(dims, f);
    const double kc_sq = (M_PI / dims.a) * (M_PI / dims.a);
    std::vector<WavePort> ports2d;
    for (int tag : {2, 3})
      ports2d.push_back(
          build_wave_port_2d(mesh, tag, mode, bc.dirichlet_edges, kc_sq));
    check_column(calculate_sparams(mesh, p, bc, ports2d, opts), kPort2d[i],
                 "calculate_sparams 2D", f);
  }

  // frequency_sweep fast path (K/M reuse), ports built at 10 GHz as in
  // the reference run.
  {
    MaxwellParams p;
    p.omega = 2 * M_PI * 10e9;
    const SweepResult sw = frequency_sweep(
        mesh, p, bc, te10_ports(mesh, dims, 10e9), kFreqs, opts);
    CHECK(sw.S_matrices.size() == kFreqs.size());
    for (size_t i = 0; i < sw.S_matrices.size(); ++i)
      check_column(sw.S_matrices[i], kTe10[i], "frequency_sweep fast",
                   kFreqs[i]);
  }

  // frequency_sweep slow path (per-frequency assembly): a PML region tag
  // that no element carries forces it without changing the matrix.
  {
    MaxwellParams p;
    p.omega = 2 * M_PI * 10e9;
    p.pml_regions.insert(-12345);
    const SweepResult sw = frequency_sweep(
        mesh, p, bc, te10_ports(mesh, dims, 10e9), kFreqs, opts);
    CHECK(sw.S_matrices.size() == kFreqs.size());
    for (size_t i = 0; i < sw.S_matrices.size(); ++i)
      check_column(sw.S_matrices[i], kTe10[i], "frequency_sweep slow",
                   kFreqs[i]);
  }

  return edgefem_test::exit_code();
}
