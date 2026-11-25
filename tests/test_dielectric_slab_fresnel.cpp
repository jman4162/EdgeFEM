// Dielectric slab at normal incidence: FEM against the analytic slab.
//
// examples/dielectric_slab.msh: a 10 x 10 x 30 mm box, eps_r = 4 slab from
// z = 13.5 to 16.5 mm (region 101), ports on z = 0 (tag 1) and z = 30 mm
// (tag 2). With PEC on the x faces (tags 5, 6) and the y faces (3, 4) left
// natural (PMC), the x-polarized TEM plane wave E = x exp(-j k z) is an
// exact mode of the box: E is normal to the PEC walls and H (along y) is
// normal to the PMC walls. The port is that mode (kc = 0), weights = the
// edge line integrals of x (the x extent of each edge), solved with the
// eigenmode port operator j beta M_s.
//
// Analytic slab (exp(+j omega t)), referred to the port planes:
//   Gamma = (1 - n) / (1 + n), P = exp(-j k0 n d),
//   r = Gamma (1 - P^2) / (1 - Gamma^2 P^2),
//   t = P (1 - Gamma^2) / (1 - Gamma^2 P^2),
//   S11 = r exp(-2 j k0 z1), S21 = t exp(-j k0 (z1 + Lz - z2)).

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
constexpr double z1 = 0.0135, z2 = 0.0165, Lz = 0.030;
constexpr double eps_slab = 4.0;

WavePort tem_port(const Mesh &mesh, int tag) {
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
    const double dx = mesh.nodes[mesh.nodeIndex.at(e.n1)].xyz.x() -
                      mesh.nodes[mesh.nodeIndex.at(e.n0)].xyz.x();
    port.weights(k) = dx;
  }
  port.mode.kc = 0.0; // TEM: beta = k0 sqrt(eps_r)
  return port;
}

void analytic(double f, std::complex<double> &S11, std::complex<double> &S21) {
  const std::complex<double> j(0.0, 1.0);
  const double k0 = 2.0 * M_PI * f / c0;
  const double n = std::sqrt(eps_slab), d = z2 - z1;
  const double g = (1.0 - n) / (1.0 + n);
  const std::complex<double> P = std::exp(-j * k0 * n * d);
  const std::complex<double> den = 1.0 - g * g * P * P;
  S11 = g * (1.0 - P * P) / den * std::exp(-2.0 * j * k0 * z1);
  S21 = P * (1.0 - g * g) / den * std::exp(-j * k0 * (z1 + Lz - z2));
}
} // namespace

int main() {
  Mesh mesh = load_gmsh_v2("examples/dielectric_slab.msh");
  BC bc = build_edge_pec(mesh, 5);                      // x = 0
  for (int e : build_edge_pec(mesh, 6).dirichlet_edges) // x = Lx
    bc.dirichlet_edges.insert(e);
  const std::vector<WavePort> ports{tem_port(mesh, 1), tem_port(mesh, 2)};

  std::cout << std::fixed << std::setprecision(4);
  for (double f : {8e9, 10e9, 12e9}) {
    MaxwellParams p;
    p.omega = 2.0 * M_PI * f;
    p.eps_r_regions[101] = eps_slab;
    const Eigen::MatrixXcd S = calculate_sparams_eigenmode(mesh, p, bc, ports);
    std::complex<double> S11_ref, S21_ref;
    analytic(f, S11_ref, S21_ref);
    const double e11 = std::abs(S(0, 0) - S11_ref);
    const double e21 = std::abs(S(1, 0) - S21_ref);
    const double ph21 = std::abs(std::arg(S(1, 0) / S21_ref)) * 180.0 / M_PI;
    std::cout << f / 1e9 << " GHz: S11 " << S(0, 0) << " (ref " << S11_ref
              << "), S21 " << S(1, 0) << " (ref " << S21_ref << "), |dS11| "
              << e11 << ", |dS21| " << e21 << ", arg err S21 " << ph21
              << " deg\n";
    // Measured on this 2 mm mesh: |dS| 0.005 / 0.020 / 0.046 and S21 phase
    // 0.4 / 1.0 / 2.6 deg at 8 / 10 / 12 GHz (dispersion error grows as
    // (k h)^2; the slab has about 6 elements per wavelength at 12 GHz). A
    // wrong sign or material gives errors near 0.5.
    CHECK(e11 < 0.06);
    CHECK(e21 < 0.06);
    CHECK(ph21 < 4.0);
    // Lossless and reciprocal.
    CHECK(std::abs(std::norm(S(0, 0)) + std::norm(S(1, 0)) - 1.0) < 0.02);
    CHECK(std::abs(S(1, 0) - S(0, 1)) < 1e-6);
  }
  return edgefem_test::exit_code();
}
