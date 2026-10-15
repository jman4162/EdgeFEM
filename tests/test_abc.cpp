// First-order ABC: the assembled term is j k0 times the boundary surface
// mass matrix, int_S N_i . N_j dS over the ABC triangles.
//
// Mesh: examples/cube_cavity.msh, the unit cube with all six faces in
// physical surface 1 and no PEC, so the ABC covers the whole boundary.
// Checks:
// - A(use_abc) - A(no abc) equals j k0 M_s entry by entry, including the
//   off-diagonal entries (the earlier diagonal-only ABC had none);
// - the sweep's frequency-independent matrix (abc_surface_matrix) is M_s;
// - Whitney edge functions reproduce a constant field exactly, so for
//   E = x-hat the quadratic form x^T M_s x is int_S |E_t|^2 dS = 4 m^2
//   (x-hat is tangential on the four faces y = 0, 1 and z = 0, 1).
#include "test_check.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

int main() {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc; // no PEC: every boundary edge carries the ABC
  constexpr double c0 = 299792458.0;
  const double k0 = 2.0;

  MaxwellParams p;
  p.omega = k0 * c0;
  const auto asmbl_noabc = assemble_maxwell(mesh, p, bc, {}, -1);
  p.use_abc = true;
  p.abc_surface_tags = {1};
  const auto asmbl_abc = assemble_maxwell(mesh, p, bc, {}, -1);

  const int m = static_cast<int>(mesh.edges.size());
  CHECK(asmbl_abc.A.rows() == m);
  const Eigen::MatrixXd Ms = Eigen::MatrixXd(
      assemble_surface_mass(mesh, p.abc_surface_tags, bc.dirichlet_edges));
  const Eigen::MatrixXcd delta =
      Eigen::MatrixXcd(asmbl_abc.A) - Eigen::MatrixXcd(asmbl_noabc.A);
  const Eigen::MatrixXcd expected =
      std::complex<double>(0.0, k0) * Ms.cast<std::complex<double>>();
  const double rel_err = (delta - expected).norm() / expected.norm();
  std::cout << "||A_abc - A_noabc - j k0 Ms|| / ||j k0 Ms|| = " << rel_err
            << std::endl;
  CHECK(rel_err < 1e-12);

  int offdiag = 0;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      if (i != j && std::abs(Ms(i, j)) > 1e-14)
        ++offdiag;
  std::cout << "off-diagonal nonzeros of Ms: " << offdiag << std::endl;
  CHECK(offdiag > 0);
  CHECK((Ms - Ms.transpose()).norm() < 1e-14 * Ms.norm());

  const Eigen::MatrixXd Ms_sweep =
      Eigen::MatrixXcd(abc_surface_matrix(mesh, p, bc)).real();
  CHECK((Ms_sweep - Ms).norm() < 1e-14 * Ms.norm());

  // Edge DOFs of the constant field E = x-hat: the line integral along the
  // global edge direction (lower node id to higher).
  Eigen::VectorXd x(m);
  for (int e = 0; e < m; ++e) {
    const auto &edge = mesh.edges[e];
    const Eigen::Vector3d r0 = mesh.nodes[mesh.nodeIndex.at(edge.n0)].xyz;
    const Eigen::Vector3d r1 = mesh.nodes[mesh.nodeIndex.at(edge.n1)].xyz;
    x(e) = r1.x() - r0.x();
  }
  const double tangential_power = x.dot(Ms * x);
  std::cout << "x^T Ms x for E = x-hat: " << tangential_power << " (expected 4)"
            << std::endl;
  CHECK_NEAR(tangential_power, 4.0, 1e-10);

  // use_abc = false: no ABC term.
  p.use_abc = false;
  CHECK(abc_surface_matrix(mesh, p, bc).nonZeros() == 0);

  return edgefem_test::exit_code();
}
