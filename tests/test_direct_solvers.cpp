// Direct solver backends, Dirichlet elimination and solve-failure contract.
//
// Modes (argv[1]):
//   (none)    elimination, failures, sweep bookkeeping; any build.
//   backends  UMFPACK and SparseLU give the same S-parameters to 1e-9
//             relative. Exits 77 (CTest: skipped) in a build without
//             SuiteSparse.
//
// Meshes: examples/parallel_plate_line.msh (lumped two-port TEM line, see
// test_lumped_tem_line.cpp) and examples/rect_waveguide.msh (WR-90 section
// with analytic TE10 wave ports, see test_frequency_sweep.cpp).
//
// Elimination: calculate_sparams_with_fields and frequency_sweep factor
// the system with the Dirichlet (PEC) rows and columns removed and expand
// the solution back to one entry per mesh edge. The reference is the
// former path: the full matrix from assemble_maxwell, whose Dirichlet rows
// and columns are identity rows, solved with the same backend. PEC entries
// must be exactly zero and the free entries equal to 1e-9 relative.

#include "edgefem/bc.hpp"
#include "edgefem/free_dofs.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/ports/lumped_port.hpp"
#include "edgefem/ports/wave_port.hpp"
#include "edgefem/solver.hpp"
#include "edgefem/sweep.hpp"
#include "test_check.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace edgefem;

namespace {

constexpr double kC0 = 299792458.0;
constexpr double kRelTol = 1e-9;
constexpr int kSkip = 77;

// ---- Cases ---------------------------------------------------------------

struct Case {
  std::string name;
  Mesh mesh;
  BC bc;
  std::vector<WavePort> ports;
  MaxwellParams p; // omega set per frequency
  std::vector<double> freqs;
};

WavePort lumped_port(const Mesh &mesh, int tag) {
  LumpedPortConfig cfg;
  cfg.surface_tag = tag;
  cfg.z0 = 50.0;
  cfg.e_direction = Eigen::Vector3d(0, 0, 1);
  return build_lumped_port(mesh, cfg);
}

// Lumped two-port on the parallel-plate line. `lossy` adds a lossy
// dielectric and the ABC on every boundary face that is not PEC.
Case tem_line(bool lossy) {
  Case c;
  c.name = lossy ? "tem_line_lossy_abc" : "tem_line";
  c.mesh = load_gmsh_v2("examples/parallel_plate_line.msh");
  c.bc = build_edge_pec(c.mesh, 1);
  c.ports = {lumped_port(c.mesh, 2), lumped_port(c.mesh, 3)};
  c.p.use_abc = lossy;
  c.p.eps_r = lossy ? std::complex<double>(2.2, -0.05) : 1.0;
  c.freqs = {2e9, 5e9, 8e9};
  return c;
}

// WR-90 section with analytic TE10 wave ports at 10 GHz.
Case waveguide() {
  Case c;
  c.name = "wr90";
  c.mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  c.bc = build_edge_pec(c.mesh, 1);
  const RectWaveguidePort dims{0.02286, 0.01016};
  const double f = 10e9;
  for (int tag : {2, 3}) {
    PortSurfaceMesh surf = extract_surface_mesh(c.mesh, tag);
    PortMode mode = solve_te10_mode(dims, f);
    populate_te10_field(surf, dims, mode);
    c.ports.push_back(build_wave_port(c.mesh, surf, mode));
  }
  c.freqs = {f};
  return c;
}

MaxwellParams at(const Case &c, double f) {
  MaxwellParams p = c.p;
  p.omega = 2.0 * M_PI * f;
  return p;
}

SolveOptions direct(DirectBackend backend) {
  SolveOptions opts;
  opts.use_direct = true;
  opts.direct_backend = backend;
  return opts;
}

double rel_diff(const Eigen::MatrixXcd &a, const Eigen::MatrixXcd &b) {
  return (a - b).norm() / b.norm();
}

double rel_diff(const VecC &a, const VecC &b) {
  return (a - b).norm() / b.norm();
}

// ---- backends -------------------------------------------------------------

int run_backends() {
  if (!umfpack_available()) {
    std::cout << "built without SuiteSparse UMFPACK: skipped" << std::endl;
    return kSkip;
  }
  const auto umf = direct(DirectBackend::UmfPack);
  const auto slu = direct(DirectBackend::SparseLU);

  for (const Case &c : {tem_line(false), tem_line(true), waveguide()}) {
    for (double f : c.freqs) {
      const MaxwellParams p = at(c, f);
      const auto a =
          calculate_sparams_with_fields(c.mesh, p, c.bc, c.ports, umf);
      const auto b =
          calculate_sparams_with_fields(c.mesh, p, c.bc, c.ports, slu);
      const double ds = rel_diff(a.S, b.S);
      std::cout << c.name << " f = " << f / 1e9
                << " GHz: |S_umfpack - S_sparselu| / |S| = " << ds << std::endl;
      CHECK(a.S.allFinite());
      CHECK(ds < kRelTol);
      for (std::size_t k = 0; k < a.x.size() && k < b.x.size(); ++k)
        CHECK(rel_diff(a.x[k], b.x[k]) < kRelTol);
    }

    const SweepResult su =
        frequency_sweep(c.mesh, c.p, c.bc, c.ports, c.freqs, umf);
    const SweepResult ss =
        frequency_sweep(c.mesh, c.p, c.bc, c.ports, c.freqs, slu);
    CHECK(su.S_matrices.size() == c.freqs.size());
    CHECK(ss.S_matrices.size() == c.freqs.size());
    for (std::size_t i = 0;
         i < su.S_matrices.size() && i < ss.S_matrices.size(); ++i) {
      const double ds = rel_diff(su.S_matrices[i], ss.S_matrices[i]);
      std::cout << c.name << " sweep point " << i
                << ": |S_umfpack - S_sparselu| / |S| = " << ds << std::endl;
      CHECK(ds < kRelTol);
    }
  }
  return edgefem_test::exit_code();
}

// ---- elimination ----------------------------------------------------------

// The former solve: full matrix with identity rows for Dirichlet edges.
std::vector<VecC> full_system_solutions(const Case &c, double f,
                                        const SolveOptions &opts) {
  const MaxwellParams p = at(c, f);
  const auto asmbl = assemble_maxwell(c.mesh, p, c.bc, c.ports, -1);
  const int m = static_cast<int>(c.mesh.edges.size());
  std::vector<VecC> rhs;
  for (const auto &port : c.ports)
    rhs.push_back(port_excitation(port, c.bc, m, p.omega / kC0));
  std::vector<VecC> xs;
  for (auto &sol : solve_linear_multi(asmbl.A, rhs, opts)) {
    CHECK(sol.converged);
    xs.push_back(sol.x);
  }
  return xs;
}

void check_same_solution(const char *what, const Case &c, const VecC &x,
                         const VecC &x_full) {
  CHECK(x.size() == x_full.size());
  if (x.size() != x_full.size())
    return;
  double pec_max = 0.0;
  for (int e : c.bc.dirichlet_edges)
    pec_max = std::max(pec_max, std::abs(x(e)));
  const double d = rel_diff(x, x_full);
  std::cout << c.name << " " << what << ": max |x_pec| = " << pec_max
            << ", |x - x_full| / |x_full| = " << d << std::endl;
  CHECK(pec_max == 0.0);
  CHECK(d < kRelTol);
}

void check_free_dofs_map() {
  // Edges 1 and 3 of 5 are Dirichlet.
  BC bc;
  bc.dirichlet_edges = {1, 3};
  const FreeDofs dofs(5, bc);
  CHECK(dofs.num_full() == 5);
  CHECK(dofs.num_free() == 3);
  CHECK(dofs.free_index(0) == 0);
  CHECK(dofs.free_index(1) == -1);
  CHECK(dofs.free_index(2) == 1);
  CHECK(dofs.free_index(4) == 2);

  Eigen::MatrixXcd dense(5, 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      dense(i, j) = std::complex<double>(10 * i + j, i - j);
  const SpMatC A = dense.sparseView();
  const Eigen::MatrixXcd Af = Eigen::MatrixXcd(dofs.restrict_matrix(A));
  CHECK(Af.rows() == 3 && Af.cols() == 3);
  const int keep[3] = {0, 2, 4};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      CHECK(Af(i, j) == dense(keep[i], keep[j]));

  VecC v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  const VecC vf = dofs.restrict_vector(v);
  CHECK(vf.size() == 3 && vf(0) == 1.0 && vf(1) == 3.0 && vf(2) == 5.0);
  const VecC back = dofs.expand(vf);
  CHECK(back.size() == 5);
  CHECK(back(0) == 1.0 && back(1) == 0.0 && back(2) == 3.0 && back(3) == 0.0 &&
        back(4) == 5.0);
}

void run_elimination() {
  check_free_dofs_map();

  for (DirectBackend backend :
       {DirectBackend::Default, DirectBackend::SparseLU}) {
    const SolveOptions opts = direct(backend);
    for (const Case &c : {tem_line(false), tem_line(true), waveguide()}) {
      const double f = c.freqs.back();
      const auto ref = full_system_solutions(c, f, opts);

      const auto sol =
          calculate_sparams_with_fields(c.mesh, at(c, f), c.bc, c.ports, opts);
      CHECK(sol.x.size() == ref.size());
      for (std::size_t k = 0; k < sol.x.size() && k < ref.size(); ++k)
        check_same_solution("calculate_sparams_with_fields", c, sol.x[k],
                            ref[k]);

      SweepFieldRequest keep;
      keep.indices = {c.freqs.size() - 1};
      const SweepResult sweep =
          frequency_sweep(c.mesh, c.p, c.bc, c.ports, c.freqs, opts, keep);
      CHECK(sweep.fields.size() == 1);
      if (sweep.fields.size() == 1) {
        const auto &xs = sweep.fields[0].x;
        CHECK(xs.size() == ref.size());
        for (std::size_t k = 0; k < xs.size() && k < ref.size(); ++k)
          check_same_solution("frequency_sweep", c, xs[k], ref[k]);
      }
    }
  }
}

// ---- failures -------------------------------------------------------------

// An iterative solve stopped after one step does not converge. The
// S-matrix column is NaN and the field is empty, as maxwell.hpp and
// sweep.hpp promise; a caller can test len(x) > 0.
void run_failures() {
  const Case c = tem_line(false);
  SolveOptions opts;
  opts.use_direct = false;
  opts.use_bicgstab = true;
  opts.use_ilut = false;
  opts.max_iterations = 1;
  opts.auto_fallback = false;

  const auto sol =
      calculate_sparams_with_fields(c.mesh, at(c, 5e9), c.bc, c.ports, opts);
  CHECK(sol.x.size() == c.ports.size());
  for (const VecC &x : sol.x)
    CHECK(x.size() == 0);
  CHECK(sol.S.array().isNaN().all());

  SweepFieldRequest keep;
  keep.indices = {0};
  keep.min_reflection = true;
  const SweepResult sweep =
      frequency_sweep(c.mesh, c.p, c.bc, c.ports, {5e9}, opts, keep);
  CHECK(sweep.fields.size() == 1);
  for (const auto &kept : sweep.fields)
    for (const VecC &x : kept.x)
      CHECK(x.size() == 0);
  // No frequency had a finite |S11|, so none is the minimum.
  CHECK(!sweep.min_reflection_index.has_value());

  // A structurally singular matrix (empty column) fails to factor with
  // every backend: converged is false and x is empty.
  for (DirectBackend backend :
       {DirectBackend::Default, DirectBackend::SparseLU}) {
    SpMatC A(3, 3);
    A.insert(0, 0) = 1.0;
    A.insert(1, 1) = 1.0;
    A.insert(2, 1) = 1.0;
    A.makeCompressed();
    const VecC b = VecC::Ones(3);
    const auto res = solve_linear_multi(A, {b}, direct(backend));
    CHECK(res.size() == 1);
    for (const auto &r : res) {
      std::cout << "singular matrix, " << r.method
                << ": converged = " << r.converged << ", error '"
                << r.error_message << "'" << std::endl;
      CHECK(!r.converged);
      CHECK(r.x.size() == 0);
      // The message carries the backend's own status.
      if (r.method == "UMFPACK")
        CHECK(r.error_message.find(
                  "status 1 (UMFPACK_WARNING_singular_matrix)") !=
              std::string::npos);
      else
        CHECK(r.error_message.find("SparseLU factorization failed: ") == 0);
    }
  }
}

// ---- sweep bookkeeping ----------------------------------------------------

void run_sweep_bookkeeping() {
  Case c = tem_line(false);
  // A 25 ohm port on the 50 ohm line reflects, so |S11| varies with
  // frequency and the minimum is a real choice.
  for (auto &port : c.ports)
    port.mode.Z0 = 25.0;
  const std::vector<double> freqs = {2e9, 3.5e9, 5e9, 6.5e9, 8e9};
  const SolveOptions opts = direct(DirectBackend::Default);

  const auto n0 = direct_factorization_count();
  const SweepResult plain =
      frequency_sweep(c.mesh, c.p, c.bc, c.ports, freqs, opts);
  CHECK(direct_factorization_count() - n0 == freqs.size());
  CHECK(!plain.min_reflection_index.has_value());
  CHECK(plain.fields.empty());

  std::size_t argmin = 0;
  for (std::size_t i = 1; i < freqs.size(); ++i)
    if (std::abs(plain.S_matrices[i](0, 0)) <
        std::abs(plain.S_matrices[argmin](0, 0)))
      argmin = i;

  SweepFieldRequest keep;
  keep.min_reflection = true;
  const SweepResult best =
      frequency_sweep(c.mesh, c.p, c.bc, c.ports, freqs, opts, keep);
  CHECK(best.min_reflection_index.has_value());
  if (best.min_reflection_index) {
    std::cout << "min |S11| at index " << *best.min_reflection_index
              << " (argmin " << argmin << ")" << std::endl;
    CHECK(*best.min_reflection_index == argmin);
  }
  CHECK(best.fields.size() == 1);
  if (best.fields.size() == 1)
    CHECK(best.fields[0].index == argmin);

  const auto n1 = direct_factorization_count();
  calculate_sparams(c.mesh, at(c, 5e9), c.bc, c.ports, opts);
  CHECK(direct_factorization_count() - n1 == 1);
}

} // namespace

int main(int argc, char **argv) {
  std::cout.precision(6);
  if (argc > 1 && std::strcmp(argv[1], "backends") == 0)
    return run_backends();
  run_elimination();
  run_failures();
  run_sweep_bookkeeping();
  return edgefem_test::exit_code();
}
