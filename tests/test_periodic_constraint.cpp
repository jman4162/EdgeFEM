// Floquet convention, edge pairing and constraint elimination of the
// periodic BC, on the 5 mm unit-cell mesh python/tests/data/
// unit_cell_bare_d8.msh (x faces tagged 5 -> 6, y faces 7 -> 8).
//
// Time convention exp(+j omega t): a field with transverse wavevector k_t
// obeys E(r + L) = exp(-j k_t . L) E(r). For a plane wave
// E = e0 exp(-j k . r) the edge DOFs (line integrals from n0 to n1) obey
// x_slave = phase_shift * orient * x_master exactly, on any mesh, when
// phase_shift and orient are right.

#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iostream>

#include "edgefem/maxwell.hpp"
#include "edgefem/periodic.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

namespace {
constexpr double kC0 = 299792458.0;
const char *kMesh = "python/tests/data/unit_cell_bare_d8.msh";

VecC plane_wave_dofs(const Mesh &mesh, const Eigen::Vector3d &k,
                     const Eigen::Vector3cd &e0) {
  const std::complex<double> j(0.0, 1.0);
  VecC x(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    const auto &p0 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[i].n0)].xyz;
    const auto &p1 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[i].n1)].xyz;
    const Eigen::Vector3d d = p1 - p0;
    const double kd = k.dot(d);
    std::complex<double> avg = std::exp(-j * k.dot(p0));
    if (std::abs(kd) > 1e-12)
      avg *= (std::exp(-j * kd) - 1.0) / (-j * kd);
    x(i) = (e0(0) * d.x() + e0(1) * d.y() + e0(2) * d.z()) * avg;
  }
  return x;
}

Eigen::Vector3d edge_dir(const Mesh &mesh, int e) {
  const auto &p0 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[e].n0)].xyz;
  const auto &p1 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[e].n1)].xyz;
  return p1 - p0;
}

struct Setup {
  Mesh mesh;
  double omega;
  Eigen::Vector3d k;
  Eigen::Vector3cd e0;
};

Setup make_setup() {
  Setup s{load_gmsh_v2(kMesh), 2 * M_PI * 20e9, {}, {}};
  const double k0 = s.omega / kC0;
  const double th = 40.0 * M_PI / 180.0, ph = 30.0 * M_PI / 180.0;
  s.k = k0 * Eigen::Vector3d(std::sin(th) * std::cos(ph),
                             std::sin(th) * std::sin(ph), std::cos(th));
  s.e0 = Eigen::Vector3cd(1.0, std::complex<double>(0.0, 0.5), -0.3);
  return s;
}

void check_bloch_relation(const Setup &s, int master, int slave,
                          const Eigen::Vector3d &L) {
  PeriodicBC pbc = build_periodic_pairs(s.mesh, master, slave, L);
  set_floquet_phase(pbc, Eigen::Vector2d(s.k.x(), s.k.y()));
  const std::complex<double> j(0.0, 1.0);
  CHECK(std::abs(pbc.phase_shift - std::exp(-j * s.k.dot(L))) < 1e-12);

  const VecC x = plane_wave_dofs(s.mesh, s.k, s.e0);
  double worst = 0.0;
  int bad_orient = 0;
  for (const auto &pr : pbc.pairs) {
    const double o = pr.master_orient * pr.slave_orient;
    const double geo =
        edge_dir(s.mesh, pr.master_edge).dot(edge_dir(s.mesh, pr.slave_edge));
    if ((geo > 0) != (o > 0))
      ++bad_orient;
    const std::complex<double> want = pbc.phase_shift * o * x(pr.master_edge);
    worst = std::max(worst, std::abs(x(pr.slave_edge) - want));
  }
  std::cout << "pairs " << master << "->" << slave << ": " << pbc.pairs.size()
            << ", wrong orientation " << bad_orient
            << ", max |x_s - phase o x_m| = " << worst << "\n";
  CHECK(!pbc.pairs.empty());
  CHECK(bad_orient == 0);
  // Gmsh writes the translated nodes to ~1e-16 relative accuracy.
  CHECK(worst < 1e-9 * x.cwiseAbs().maxCoeff());
}

// Eliminating the slave edges and expanding the solution reproduces a
// Bloch field exactly: with b = A x_true for a field that satisfies the
// constraint, the reduced solve must return x_true.
void check_reduce_expand_round_trip(const Setup &s) {
  const Eigen::Vector3d Lx(5e-3, 0.0, 0.0);
  PeriodicBC pbc = build_periodic_pairs(s.mesh, 5, 6, Lx);
  set_floquet_phase(pbc, Eigen::Vector2d(s.k.x(), s.k.y()));

  BC bc; // no PEC
  MaxwellParams p;
  p.omega = s.omega;
  const VecC x_true = plane_wave_dofs(s.mesh, s.k, s.e0);
  MaxwellAssembly sys = assemble_maxwell(s.mesh, p, bc, {}, -1);
  VecC b = sys.A * x_true;
  eliminate_periodic(sys.A, b, pbc, bc);

  SolveOptions opts;
  opts.use_direct = true;
  SolveResult res = solve_linear(sys.A, b, opts);
  CHECK(res.converged);
  const VecC x = expand_periodic(res.x, pbc, bc);
  const double err = (x - x_true).norm() / x_true.norm();
  std::cout << "reduce/expand round trip: relative error " << err << "\n";
  CHECK(err < 1e-8);
}
} // namespace

int main() {
  const Setup s = make_setup();
  check_bloch_relation(s, 5, 6, Eigen::Vector3d(5e-3, 0.0, 0.0));
  check_bloch_relation(s, 7, 8, Eigen::Vector3d(0.0, 5e-3, 0.0));
  check_reduce_expand_round_trip(s);
  return edgefem_test::exit_code();
}
