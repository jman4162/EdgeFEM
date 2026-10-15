#include "test_check.hpp"
#include <Eigen/Dense>
#include <string>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

void test_symmetry() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);
  CHECK(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  Eigen::MatrixXcd A = Eigen::MatrixXcd(asmbl.A);
  CHECK(A.isApprox(A.adjoint(), 1e-9));
}

// Gap-port source and load (derivation in src/assemble_maxwell.cpp) for a
// Lumped port. The port's weights are its voltage functional v (here v = 1
// on one edge).
// With omega = 1 rad/s, k0 eta0 = omega mu0 = 4 pi 1e-7 exactly, so
//   source  b = j k0 eta0 (2 / sqrt(Z0)) v = j 3.5543e-7  (Z0 = 50 ohm),
//   load    A += j k0 eta0 / Z0 v v^H        = j 2.5133e-8 on the diagonal.
// Until 2026-10 this test expected b = 2 / sqrt(Z0) = 0.2828, the source of
// the non-physical port term w w^H / Z0 (a real load).
void test_wave_port_excitation() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;

  WavePort port1;
  port1.edges = {10};
  port1.weights = Eigen::VectorXcd::Ones(1);
  port1.mode.Z0 = 50.0;
  port1.kind = PortKind::Lumped;

  WavePort port2 = port1;
  port2.edges[0] = 20;

  // The checks below need port edges that are not PEC.
  CHECK(!bc.dirichlet_edges.count(port1.edges[0]));
  CHECK(!bc.dirichlet_edges.count(port2.edges[0]));

  // Excite port 1
  std::vector<WavePort> ports{port1, port2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);
  auto bare = assemble_maxwell(mesh, p, bc, {}, -1);

  const double Z0 = 50.0;
  const double mu0 = 4.0 * M_PI * 1e-7;
  const double k0_eta0 = p.omega * mu0;
  const std::complex<double> j(0.0, 1.0);
  const std::complex<double> b_ref = j * k0_eta0 * 2.0 / std::sqrt(Z0);
  const std::complex<double> y_ref = j * k0_eta0 / Z0;

  // Source on the active port's edge only.
  CHECK(std::abs(asmbl.b(port1.edges[0]) - b_ref) < 1e-9 * std::abs(b_ref));
  CHECK(std::abs(asmbl.b(port2.edges[0])) == 0.0);
  CHECK(std::abs(asmbl.b.norm() - std::abs(b_ref)) < 1e-9 * std::abs(b_ref));

  // Load j k0 eta0 / Z0 on each port's diagonal, and nowhere else.
  for (int e : {port1.edges[0], port2.edges[0]}) {
    const std::complex<double> dy = asmbl.A.coeff(e, e) - bare.A.coeff(e, e);
    CHECK(std::abs(dy - y_ref) < 1e-9 * std::abs(y_ref));
  }
  const SpMatC diff = asmbl.A - bare.A;
  int changed = 0;
  for (int col = 0; col < diff.outerSize(); ++col)
    for (SpMatC::InnerIterator it(diff, col); it; ++it)
      if (std::abs(it.value()) > 0.0)
        ++changed;
  CHECK(changed == 2);
}

// A Modal port (the default kind; build_wave_port and the other wave-port
// builders) keeps the 1.0.0 terms: source b = 2 w / sqrt(Z0) = 0.28284 and
// load A += w w^H / Z0 = 0.02 on the diagonal, both real and without the
// j k0 eta0 factor of the gap port.
void test_modal_port_terms() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 2 * M_PI * 1e9;

  WavePort port1;
  port1.edges = {10};
  port1.weights = Eigen::VectorXcd::Ones(1);
  port1.mode.Z0 = 50.0;
  CHECK(port1.kind == PortKind::Modal);
  WavePort port2 = port1;
  port2.edges[0] = 20;
  CHECK(!bc.dirichlet_edges.count(port1.edges[0]));
  CHECK(!bc.dirichlet_edges.count(port2.edges[0]));

  std::vector<WavePort> ports{port1, port2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);
  auto bare = assemble_maxwell(mesh, p, bc, {}, -1);

  const std::complex<double> b_ref = 2.0 / std::sqrt(50.0);
  const std::complex<double> y_ref = 1.0 / 50.0;
  CHECK(std::abs(asmbl.b(port1.edges[0]) - b_ref) < 1e-12);
  CHECK(std::abs(asmbl.b(port2.edges[0])) == 0.0);
  CHECK(std::abs(asmbl.b.norm() - std::abs(b_ref)) < 1e-12);
  for (int e : {port1.edges[0], port2.edges[0]}) {
    const std::complex<double> dy = asmbl.A.coeff(e, e) - bare.A.coeff(e, e);
    CHECK(std::abs(dy - y_ref) < 1e-12);
  }
  const SpMatC diff = asmbl.A - bare.A;
  int changed = 0;
  for (int col = 0; col < diff.outerSize(); ++col)
    for (SpMatC::InnerIterator it(diff, col); it; ++it)
      if (std::abs(it.value()) > 0.0)
        ++changed;
  CHECK(changed == 2);
}

void test_sparam_calculation() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 2 * M_PI * 1e9; // 1 GHz

  WavePort port1;
  port1.edges = {10};
  port1.weights = Eigen::VectorXcd::Ones(1);
  port1.mode.Z0 = 50.0;
  port1.kind = PortKind::Lumped;

  WavePort port2 = port1;
  port2.edges[0] = 20;
  std::vector<WavePort> ports{port1, port2};

  auto S = calculate_sparams(mesh, p, bc, ports);

  CHECK(S.rows() == 2);
  CHECK(S.cols() == 2);

  // Check for passivity (for the first column)
  double s11_mag_sq = std::norm(S(0, 0));
  double s21_mag_sq = std::norm(S(1, 0));
  CHECK(s11_mag_sq + s21_mag_sq <=
        1.0 + 1e-6); // Add tolerance for numerical error
}

void test_mixed_bc() {
  Mesh mixed_mesh = load_gmsh_v2("examples/mixed_cavity.msh");
  BC mixed_bc = build_edge_pec(mixed_mesh, 1); // PEC tag is 1

  Mesh pec_mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC pec_bc = build_edge_pec(pec_mesh, 1);

  // The mixed cavity has one PMC wall, so it should have fewer PEC edges
  // than the fully PEC cavity.
  CHECK(mixed_bc.dirichlet_edges.size() < pec_bc.dirichlet_edges.size());
  CHECK(mixed_bc.dirichlet_edges.size() > 0);
}

// A lossy dielectric (eps_r = 1 - 0.1 j) adds +j 0.1 k0^2 M to A, which is
// anti-Hermitian, so A is no longer Hermitian. At 10 GHz, k0^2 = 4.4e4
// m^-2. Until 2026-10 the test used omega = 1 rad/s, where 0.1 k0^2 is
// 1.1e-18 and the term vanished below the 1e-9 isApprox tolerance.
void test_lossy_material() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 2 * M_PI * 10e9;

  MaxwellParams lossless = p;
  lossless.eps_r = {1.0, 0.0};
  p.eps_r = {1.0, -0.1}; // Lossy dielectric

  const Eigen::MatrixXcd A0 =
      Eigen::MatrixXcd(assemble_maxwell(mesh, lossless, bc, {}, -1).A);
  const Eigen::MatrixXcd A =
      Eigen::MatrixXcd(assemble_maxwell(mesh, p, bc, {}, -1).A);
  // The lossless matrix is Hermitian; the lossy one is not.
  CHECK(A0.isApprox(A0.adjoint(), 1e-9));
  CHECK(!A.isApprox(A.adjoint(), 1e-9));
  // The difference is the loss term alone: purely imaginary, with
  // Im(A - A0) = 0.1 k0^2 M >= 0 on the diagonal (zero on PEC rows).
  const Eigen::MatrixXcd dA = A - A0;
  CHECK(dA.real().cwiseAbs().maxCoeff() < 1e-9 * A0.cwiseAbs().maxCoeff());
  CHECK(dA.imag().diagonal().minCoeff() >= 0.0);
  CHECK(dA.imag().diagonal().maxCoeff() > 0.0);
}

int main() {
  test_symmetry();
  test_wave_port_excitation();
  test_modal_port_terms();
  test_sparam_calculation();
  test_mixed_bc();
  test_lossy_material();
  return edgefem_test::exit_code();
}
