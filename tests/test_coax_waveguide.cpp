// Air coaxial line (examples/coax_50ohm.msh: a = 1 mm, b = 2.3 mm,
// L = 30 mm; PEC tags 1 (inner) and 2 (outer), ports 3 at z = 0 and 4 at
// z = L) with a TEM
// wave port, through the eigenmode port operator j beta M_s (kc = 0).
//
// The TEM field is E = r_hat / r (any amplitude), so the port weight of an
// edge from p0 to p1 is its exact line integral int dr / r = ln(r1 / r0),
// with r the distance from the axis. A matched air line gives |S11| = 0,
// |S21| = 1 and phase(S21) = -k0 L at every frequency (no cutoff).

#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <set>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

namespace {
constexpr double c0 = 299792458.0;
constexpr double L = 0.030;

WavePort coax_tem_port(const Mesh &mesh, int tag) {
  std::set<int> edge_set;
  for (const auto &tri : mesh.tris)
    if (tri.phys == tag)
      for (int e : tri.edges)
        edge_set.insert(e);
  WavePort port;
  port.surface_tag = tag;
  port.edges.assign(edge_set.begin(), edge_set.end());
  port.weights.resize(port.edges.size());
  for (size_t k = 0; k < port.edges.size(); ++k) {
    const Edge &e = mesh.edges[port.edges[k]];
    const auto &p0 = mesh.nodes[mesh.nodeIndex.at(e.n0)].xyz;
    const auto &p1 = mesh.nodes[mesh.nodeIndex.at(e.n1)].xyz;
    port.weights(k) =
        std::log(std::hypot(p1.x(), p1.y()) / std::hypot(p0.x(), p0.y()));
  }
  port.mode.kc = 0.0;
  return port;
}
} // namespace

int main() {
  Mesh mesh = load_gmsh_v2("examples/coax_50ohm.msh");
  BC bc = build_edge_pec(mesh, 1);
  for (int e : build_edge_pec(mesh, 2).dirichlet_edges)
    bc.dirichlet_edges.insert(e);
  const std::vector<WavePort> ports{coax_tem_port(mesh, 3),
                                    coax_tem_port(mesh, 4)};

  std::cout << std::fixed << std::setprecision(4);
  for (double f : {2e9, 5e9, 10e9}) {
    MaxwellParams p;
    p.omega = 2.0 * M_PI * f;
    const Eigen::MatrixXcd S = calculate_sparams_eigenmode(mesh, p, bc, ports);
    const double k0 = p.omega / c0;
    const double ph_err =
        std::abs(
            std::arg(S(1, 0) * std::exp(std::complex<double>(0, k0 * L)))) *
        180.0 / M_PI;
    std::cout << f / 1e9 << " GHz: |S11| " << std::abs(S(0, 0)) << ", |S21| "
              << std::abs(S(1, 0)) << ", phase(S21) error " << ph_err
              << " deg\n";
    // Measured: |S11| 0.028 / 0.005 / 0.006, |S21| >= 0.9995, phase error
    // 0.7 / 1.8 / 3.9 deg at 2 / 5 / 10 GHz (mesh dispersion over one
    // wavelength of line at 10 GHz).
    CHECK(std::abs(S(0, 0)) < 0.05);
    CHECK(std::abs(std::abs(S(1, 0)) - 1.0) < 0.02);
    CHECK(ph_err < 6.0);
    CHECK(std::abs(S(1, 0) - S(0, 1)) < 1e-6);
  }
  return edgefem_test::exit_code();
}
