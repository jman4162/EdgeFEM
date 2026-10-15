// Factor + solve timings for the direct solver backends and orderings.
//
//   bench_direct_lu <mesh.msh> <case> <freq_hz> <config>
//
// case:   patch  StackedPatchDesign mesh (PEC tags 1, 4, 10; lumped port on
//                tag 5, E along x, 50 ohm; ABC on tag 50; substrate tag 110
//                with eps_r 3.66 (1 - j 0.0037), the RO4350B patch)
//         wr90   examples/rect_waveguide.msh with TE10 wave ports (tags 2, 3)
// config: <backend>[-full]
//         backend: sparselu-colamd, sparselu-amd, sparselu-metis (only in a
//                  build that found METIS), umfpack (UMFPACK default
//                  ordering), umfpack-metis (METIS), umfpack-guard
//                  (METIS_GUARD, what DirectLU uses); umfpack* need
//                  SuiteSparse.
//         -full:   factor the full matrix with identity rows for PEC edges
//                  (the path before Dirichlet elimination) instead of the
//                  reduced one.
//
// Prints one line: case, config, n, nnz, analyze/factor/solve seconds,
// S11, and the peak resident set size of the process. Run one config per
// process so the peak memory belongs to that config. Built only with
// EDGEFEM_BUILD_DEBUG_TOOLS; not a CTest test.

#include "edgefem/bc.hpp"
#include "edgefem/free_dofs.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/ports/lumped_port.hpp"
#include "edgefem/ports/wave_port.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseLU>
#include <iostream> // before MetisSupport.h, which uses std::cerr
#ifdef EDGEFEM_BENCH_METIS
#include <Eigen/MetisSupport>
#endif
#ifdef EDGEFEM_BENCH_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

#include <sys/resource.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace edgefem;

namespace {

constexpr double kC0 = 299792458.0;

double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

double peak_rss_gb() {
  rusage ru{};
  getrusage(RUSAGE_SELF, &ru);
  return static_cast<double>(ru.ru_maxrss) / 1e9; // bytes on macOS
}

struct Timing {
  double analyze = 0, factor = 0, solve = 0;
  bool ok = false;
  VecC x;
};

template <typename Solver>
Timing run(Solver &solver, const SpMatC &A, const VecC &b) {
  Timing t;
  double t0 = now();
  solver.analyzePattern(A);
  double t1 = now();
  solver.factorize(A);
  double t2 = now();
  if (solver.info() != Eigen::Success)
    return t;
  t.x = solver.solve(b);
  double t3 = now();
  t.analyze = t1 - t0;
  t.factor = t2 - t1;
  t.solve = t3 - t2;
  t.ok = solver.info() == Eigen::Success;
  return t;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: bench_direct_lu <mesh.msh> <patch|wr90> <freq_hz> "
                 "<config>\n";
    return 2;
  }
  const std::string mesh_path = argv[1], which = argv[2], config = argv[4];
  const double f = std::stod(argv[3]);
  const bool full =
      config.size() > 5 && config.compare(config.size() - 5, 5, "-full") == 0;
  const std::string backend =
      full ? config.substr(0, config.size() - 5) : config;

  Mesh mesh = load_gmsh_v2(mesh_path);
  BC bc;
  std::vector<WavePort> ports;
  MaxwellParams p;
  p.omega = 2.0 * M_PI * f;
  if (which == "patch") {
    for (int tag : {1, 4, 10}) {
      BC t = build_edge_pec(mesh, tag);
      bc.dirichlet_edges.insert(t.dirichlet_edges.begin(),
                                t.dirichlet_edges.end());
    }
    LumpedPortConfig cfg;
    cfg.surface_tag = 5;
    cfg.z0 = 50.0;
    cfg.e_direction = Eigen::Vector3d(1, 0, 0);
    ports.push_back(build_lumped_port(mesh, cfg));
    p.use_abc = true;
    p.abc_surface_tags = {50};
    p.eps_r_regions[110] = std::complex<double>(3.66, -3.66 * 0.0037);
  } else if (which == "wr90") {
    bc = build_edge_pec(mesh, 1);
    const RectWaveguidePort dims{0.02286, 0.01016};
    for (int tag : {2, 3}) {
      PortSurfaceMesh surf = extract_surface_mesh(mesh, tag);
      PortMode mode = solve_te10_mode(dims, f);
      populate_te10_field(surf, dims, mode);
      ports.push_back(build_wave_port(mesh, surf, mode));
    }
  } else {
    std::cerr << "unknown case " << which << "\n";
    return 2;
  }

  const int m = static_cast<int>(mesh.edges.size());
  const auto asmbl = assemble_maxwell(mesh, p, bc, ports, -1);
  const VecC b_full = port_excitation(ports[0], bc, m, p.omega / kC0);
  const FreeDofs dofs(m, bc);
  const SpMatC A = full ? asmbl.A : dofs.restrict_matrix(asmbl.A);
  const VecC b = full ? b_full : dofs.restrict_vector(b_full);
  const double rss_before = peak_rss_gb();

  Timing t;
  if (backend == "sparselu-colamd") {
    Eigen::SparseLU<SpMatC, Eigen::COLAMDOrdering<int>> s;
    t = run(s, A, b);
  } else if (backend == "sparselu-amd") {
    Eigen::SparseLU<SpMatC, Eigen::AMDOrdering<int>> s;
    t = run(s, A, b);
#ifdef EDGEFEM_BENCH_METIS
  } else if (backend == "sparselu-metis") {
    Eigen::SparseLU<SpMatC, Eigen::MetisOrdering<int>> s;
    t = run(s, A, b);
#endif
#ifdef EDGEFEM_BENCH_UMFPACK
  } else if (backend == "umfpack" || backend == "umfpack-metis" ||
             backend == "umfpack-guard") {
    Eigen::UmfPackLU<SpMatC> s;
    if (backend == "umfpack-metis")
      s.umfpackControl()(UMFPACK_ORDERING) = UMFPACK_ORDERING_METIS;
    if (backend == "umfpack-guard")
      s.umfpackControl()(UMFPACK_ORDERING) = UMFPACK_ORDERING_METIS_GUARD;
    t = run(s, A, b);
#endif
  } else {
    std::cerr << "config " << config << " is not available in this build\n";
    return 2;
  }
  if (!t.ok) {
    std::cerr << config << ": factorization or solve failed\n";
    return 1;
  }

  const VecC x = full ? t.x : dofs.expand(t.x);
  Eigen::MatrixXcd S(ports.size(), 1);
  fill_sparam_column(S, 0, ports, bc, x);
  std::printf("%s %s n=%ld nnz=%ld analyze=%.3f factor=%.3f solve=%.3f "
              "S11=%.9f%+.9fj rss_before_gb=%.3f peak_rss_gb=%.3f\n",
              which.c_str(), config.c_str(), static_cast<long>(A.rows()),
              static_cast<long>(A.nonZeros()), t.analyze, t.factor, t.solve,
              S(0, 0).real(), S(0, 0).imag(), rss_before, peak_rss_gb());
  return 0;
}
