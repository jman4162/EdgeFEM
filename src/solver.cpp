#include "edgefem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#ifdef EDGEFEM_HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#endif
#include <atomic>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace edgefem {

namespace {
std::atomic<std::uint64_t> g_direct_factorizations{0};

#ifdef EDGEFEM_HAVE_UMFPACK
constexpr bool kHaveUmfPack = true;
#else
constexpr bool kHaveUmfPack = false;
#endif

DirectBackend resolve(DirectBackend backend) {
  if (backend == DirectBackend::Default)
    return kHaveUmfPack ? DirectBackend::UmfPack : DirectBackend::SparseLU;
  if (backend == DirectBackend::UmfPack && !kHaveUmfPack)
    throw std::invalid_argument(
        "DirectBackend::UmfPack: EdgeFEM was built without SuiteSparse "
        "UMFPACK (configure with -DEDGEFEM_USE_SUITESPARSE=ON)");
  return backend;
}
} // namespace

bool umfpack_available() { return kHaveUmfPack; }

std::uint64_t direct_factorization_count() {
  return g_direct_factorizations.load();
}

// ---- DirectLU --------------------------------------------------------------

struct DirectLU::Impl {
  DirectBackend backend;
  std::string method;
  bool analyzed = false;
  bool factored = false;
  // COLAMD: on the 10 GHz patch (25k free edges) Eigen's SparseLU factors
  // in 2.2 s with it and in 478 s with AMDOrdering (tests/bench_direct_lu).
  Eigen::SparseLU<SpMatC, Eigen::COLAMDOrdering<int>> sparse_lu;
#ifdef EDGEFEM_HAVE_UMFPACK
  // UmfPackLU keeps a reference to the matrix it factored and reads it again
  // in solve() (iterative refinement), so it factors this copy.
  SpMatC umf_matrix;
  Eigen::UmfPackLU<SpMatC> umfpack;
#endif
};

DirectLU::DirectLU(DirectBackend backend) : impl_(std::make_unique<Impl>()) {
  impl_->backend = resolve(backend);
  impl_->method =
      impl_->backend == DirectBackend::UmfPack ? "UMFPACK" : "SparseLU";
}

DirectLU::~DirectLU() = default;

const std::string &DirectLU::method() const { return impl_->method; }

void DirectLU::analyze_pattern(const SpMatC &A) {
  impl_->factored = false;
#ifdef EDGEFEM_HAVE_UMFPACK
  if (impl_->backend == DirectBackend::UmfPack) {
    impl_->umf_matrix = A;
    impl_->umf_matrix.makeCompressed();
    // METIS ordering: on the 10 GHz patch at density 12 (58k free edges)
    // UMFPACK factors in 0.74 s and 1.0 GB with it, 1.80 s and 2.0 GB with
    // its default (AMD); analysis is 0.35 s instead of 0.08 s and runs once
    // per sweep. METIS_GUARD is METIS for structurally symmetric matrices
    // and COLAMD when the unsymmetric strategy meets a dense row. If this
    // UMFPACK has no METIS, analyze again with AMD.
    impl_->umfpack.umfpackControl()(UMFPACK_ORDERING) =
        UMFPACK_ORDERING_METIS_GUARD;
    impl_->umfpack.analyzePattern(impl_->umf_matrix);
    if (impl_->umfpack.info() != Eigen::Success) {
      impl_->umfpack.umfpackControl()(UMFPACK_ORDERING) = UMFPACK_ORDERING_AMD;
      impl_->umfpack.analyzePattern(impl_->umf_matrix);
    }
    impl_->analyzed = impl_->umfpack.info() == Eigen::Success;
    return;
  }
#endif
  impl_->sparse_lu.analyzePattern(A);
  impl_->analyzed = true;
}

bool DirectLU::factorize(const SpMatC &A) {
  if (!impl_->analyzed)
    analyze_pattern(A);
  ++g_direct_factorizations;
  impl_->factored = false;
  if (!impl_->analyzed)
    return false;
#ifdef EDGEFEM_HAVE_UMFPACK
  if (impl_->backend == DirectBackend::UmfPack) {
    impl_->umf_matrix = A;
    impl_->umf_matrix.makeCompressed();
    impl_->umfpack.factorize(impl_->umf_matrix);
    impl_->factored = impl_->umfpack.info() == Eigen::Success;
    return impl_->factored;
  }
#endif
  impl_->sparse_lu.factorize(A);
  impl_->factored = impl_->sparse_lu.info() == Eigen::Success;
  return impl_->factored;
}

bool DirectLU::solve(const VecC &b, VecC &x) const {
  x.resize(0);
  if (!impl_->factored)
    return false;
  VecC y;
  bool ok = false;
#ifdef EDGEFEM_HAVE_UMFPACK
  if (impl_->backend == DirectBackend::UmfPack) {
    y = impl_->umfpack.solve(b);
    ok = impl_->umfpack.info() == Eigen::Success;
  } else
#endif
  {
    y = impl_->sparse_lu.solve(b);
    ok = impl_->sparse_lu.info() == Eigen::Success;
  }
  if (ok)
    x = std::move(y);
  return ok;
}

std::vector<SolveResult>
DirectLU::factorize_and_solve(const SpMatC &A, const std::vector<VecC> &rhs) {
  const bool factored = factorize(A);
  std::vector<SolveResult> out;
  out.reserve(rhs.size());
  for (const VecC &b : rhs) {
    SolveResult res;
    res.method = method();
    if (!factored) {
      res.error_message = method() + " factorization failed";
    } else if (!solve(b, res.x)) {
      res.error_message = method() + " solve failed";
    } else {
      res.converged = true;
      res.iters = 1;
      const double b_norm = b.norm();
      res.residual = b_norm > 0.0 ? (A * res.x - b).norm() / b_norm : 0.0;
    }
    out.push_back(std::move(res));
  }
  return out;
}

// ---- Direct solves ---------------------------------------------------------

namespace {

// One factorization, one solve per right-hand side.
std::vector<SolveResult> solve_direct_multi(const SpMatC &A,
                                            const std::vector<VecC> &rhs,
                                            DirectBackend backend,
                                            bool verbose) {
  DirectLU lu(backend);
  auto out = lu.factorize_and_solve(A, rhs);
  if (verbose && !out.empty() && !out.front().converged)
    std::cerr << "Solver: " << out.front().error_message << "\n";
  return out;
}

SolveResult solve_direct(const SpMatC &A, const VecC &b, DirectBackend backend,
                         bool verbose) {
  return std::move(solve_direct_multi(A, {b}, backend, verbose).front());
}

} // namespace

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt) {
  SolveResult res;

  // Print header if verbose
  if (opt.verbose) {
    const char *method = opt.use_direct     ? "direct LU"
                         : opt.use_bicgstab ? "BiCGSTAB"
                                            : "CG";
    const char *precond = (!opt.use_direct && opt.use_ilut) ? "+ILUT" : "";
    std::cerr << "Solver: Starting " << method << precond
              << " solve (N=" << A.rows() << ", nnz=" << A.nonZeros() << ")\n";
    std::cerr << std::fixed << std::setprecision(6);
  }

  if (opt.use_direct) {
    res = solve_direct(A, b, opt.direct_backend, opt.verbose);
  } else if (opt.use_bicgstab) {
    if (opt.use_ilut) {
      // BiCGSTAB with ILUT preconditioner
      Eigen::BiCGSTAB<SpMatC, Eigen::IncompleteLUT<std::complex<double>>>
          solver;
      solver.preconditioner().setFillfactor(
          static_cast<int>(opt.ilut_fill_factor));
      solver.preconditioner().setDroptol(opt.ilut_drop_tolerance);
      solver.setTolerance(opt.tolerance);
      solver.setMaxIterations(opt.max_iterations);
      solver.compute(A);
      if (solver.info() != Eigen::Success) {
        res.method = "BiCGSTAB+ILUT";
        res.converged = false;
        res.error_message = "ILUT preconditioner setup failed";
        if (opt.verbose) {
          std::cerr << "Solver: ILUT preconditioner setup FAILED\n";
        }
        // Fall through to auto-fallback below
      } else {
        res.x = solver.solve(b);
        res.method = "BiCGSTAB+ILUT";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
        res.converged = (solver.info() == Eigen::Success);
        if (!res.converged) {
          res.error_message = "Solver did not converge within max iterations";
        }
      }
    } else {
      // BiCGSTAB with diagonal (default) preconditioner
      Eigen::BiCGSTAB<SpMatC> solver;
      solver.setTolerance(opt.tolerance);
      solver.setMaxIterations(opt.max_iterations);
      solver.compute(A);
      if (solver.info() != Eigen::Success) {
        res.method = "BiCGSTAB";
        res.converged = false;
        res.error_message = "Matrix decomposition failed";
        if (opt.verbose) {
          std::cerr << "Solver: Matrix decomposition FAILED\n";
        }
      } else {
        res.x = solver.solve(b);
        res.method = "BiCGSTAB";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
        res.converged = (solver.info() == Eigen::Success);
        if (!res.converged) {
          res.error_message = "Solver did not converge within max iterations";
        }
      }
    }

    // Auto-fallback to direct solver if iterative solver failed
    if (!res.converged && opt.auto_fallback) {
      std::string failed_method = res.method;
      if (opt.verbose) {
        std::cerr << "Solver: " << failed_method
                  << " failed: falling back to a direct solver\n";
      }
      res = solve_direct(A, b, opt.direct_backend, opt.verbose);
      res.method = failed_method + "->" + res.method;
    }
  } else {
    // CG solver
    if (opt.use_ilut) {
      Eigen::ConjugateGradient<SpMatC, Eigen::Lower | Eigen::Upper,
                               Eigen::IncompleteLUT<std::complex<double>>>
          solver;
      solver.preconditioner().setFillfactor(
          static_cast<int>(opt.ilut_fill_factor));
      solver.preconditioner().setDroptol(opt.ilut_drop_tolerance);
      solver.setTolerance(opt.tolerance);
      solver.setMaxIterations(opt.max_iterations);
      solver.compute(A);
      if (solver.info() != Eigen::Success) {
        res.method = "CG+ILUT";
        res.converged = false;
        res.error_message = "ILUT preconditioner setup failed";
      } else {
        res.x = solver.solve(b);
        res.method = "CG+ILUT";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
        res.converged = (solver.info() == Eigen::Success);
        if (!res.converged) {
          res.error_message = "Solver did not converge within max iterations";
        }
      }
    } else {
      Eigen::ConjugateGradient<SpMatC> solver;
      solver.setTolerance(opt.tolerance);
      solver.setMaxIterations(opt.max_iterations);
      solver.compute(A);
      if (solver.info() != Eigen::Success) {
        res.method = "CG";
        res.converged = false;
        res.error_message = "Matrix decomposition failed";
        if (opt.verbose) {
          std::cerr << "Solver: Matrix decomposition FAILED\n";
        }
      } else {
        res.x = solver.solve(b);
        res.method = "CG";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
        res.converged = (solver.info() == Eigen::Success);
        if (!res.converged) {
          res.error_message = "Solver did not converge within max iterations";
        }
      }
    }

    // Auto-fallback for CG too
    if (!res.converged && opt.auto_fallback) {
      std::string failed_method = res.method;
      if (opt.verbose) {
        std::cerr << "Solver: " << failed_method
                  << " failed: falling back to a direct solver\n";
      }
      res = solve_direct(A, b, opt.direct_backend, opt.verbose);
      res.method = failed_method + "->" + res.method;
    }
  }

  // Print summary if verbose
  if (opt.verbose) {
    std::cerr << "Solver: " << res.method << " "
              << (res.converged ? "CONVERGED" : "FAILED") << " in " << res.iters
              << " iterations"
              << ", residual=" << std::scientific << res.residual << "\n";
    std::cerr << std::fixed;
  }

  // Call progress callback with final state if provided
  if (opt.progress_callback) {
    opt.progress_callback(res.iters, res.residual);
  }

  return res;
}

std::vector<SolveResult> solve_linear_multi(const SpMatC &A,
                                            const std::vector<VecC> &rhs,
                                            const SolveOptions &opt) {
  if (opt.use_direct)
    return solve_direct_multi(A, rhs, opt.direct_backend, opt.verbose);
  std::vector<SolveResult> out;
  out.reserve(rhs.size());
  for (const VecC &b : rhs)
    out.push_back(solve_linear(A, b, opt));
  return out;
}

} // namespace edgefem
