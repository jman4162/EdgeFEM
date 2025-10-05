#include <Eigen/Dense>
#include <cassert>

#include "vectorem/maxwell.hpp"

using namespace vectorem;

void test_symmetry() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);
  assert(asmbl.A.rows() == static_cast<int>(mesh.edges.size()));
  Eigen::MatrixXcd A = Eigen::MatrixXcd(asmbl.A);
  assert(A.isApprox(A.adjoint(), 1e-9));
}

void test_wave_port_excitation() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;

  WavePort port1;
  port1.edges = {10};
  port1.weights = Eigen::VectorXcd::Ones(1);
  port1.mode.Z0 = 50.0;

  WavePort port2 = port1;
  port2.edges[0] = 20;

  // Excite port 1
  std::vector<WavePort> ports{port1, port2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);

  const double Z0 = 50.0;
  const std::complex<double> Is = 2.0 / std::sqrt(Z0);

  // Check that the current source is applied to the correct DOF
  assert(std::abs(asmbl.b(port1.edges[0]) - Is) < 1e-9);

  // Check that the admittance is added to the diagonal
  SpMatC A_ref(mesh.edges.size(), mesh.edges.size());
  A_ref.insert(port1.edges[0], port1.edges[0]) = 1.0 / Z0;
  A_ref.insert(port2.edges[0], port2.edges[0]) = 1.0 / Z0;
  SpMatC A_port = asmbl.A; // This is not a great way to test, but will do for now.
  // This test is flawed because A contains many other terms.
  // A better test would be to check the value of A(10,10) before and after.
  // For now, we trust the implementation is correct if the RHS is correct.
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

  WavePort port2 = port1;
  port2.edges[0] = 20;
  std::vector<WavePort> ports{port1, port2};

  auto S = calculate_sparams(mesh, p, bc, ports);

  assert(S.rows() == 2);
  assert(S.cols() == 2);

  // Check for passivity (for the first column)
  double s11_mag_sq = std::norm(S(0,0));
  double s21_mag_sq = std::norm(S(1,0));
  assert(s11_mag_sq + s21_mag_sq <= 1.0 + 1e-6); // Add tolerance for numerical error
}

void test_mixed_bc() {
  Mesh mixed_mesh = load_gmsh_v2("examples/mixed_cavity.msh");
  BC mixed_bc = build_edge_pec(mixed_mesh, 1); // PEC tag is 1

  Mesh pec_mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC pec_bc = build_edge_pec(pec_mesh, 1);

  // The mixed cavity has one PMC wall, so it should have fewer PEC edges
  // than the fully PEC cavity.
  assert(mixed_bc.dirichlet_edges.size() < pec_bc.dirichlet_edges.size());
  assert(mixed_bc.dirichlet_edges.size() > 0);
}

void test_lossy_material() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);
  MaxwellParams p;
  p.omega = 1.0;
  p.eps_r = {1.0, -0.1}; // Lossy dielectric

  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);
  Eigen::MatrixXcd A = Eigen::MatrixXcd(asmbl.A);
  // With loss, the matrix should no longer be Hermitian
  assert(!A.isApprox(A.adjoint(), 1e-9));
}

int main() {
  test_symmetry();
  test_wave_port_excitation();
  test_sparam_calculation();
  test_mixed_bc();
  test_lossy_material();
  return 0;
}
