#include "edgefem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <iomanip>
#include <iostream>

namespace edgefem {

// SparseLU direct solve (shared between direct path and fallback)
static SolveResult solve_direct(const SpMatC &A, const VecC &b, bool verbose) {
  SolveResult res;
  res.method = "SparseLU";

  Eigen::SparseLU<SpMatC> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    res.converged = false;
    res.error_message = "SparseLU factorization failed";
    if (verbose) {
      std::cerr << "Solver: SparseLU factorization FAILED\n";
    }
    return res;
  }
  res.x = solver.solve(b);
  res.iters = 1;
  res.residual = (A * res.x - b).norm() / b.norm();
  res.converged = (solver.info() == Eigen::Success);
  if (!res.converged) {
    res.error_message = "SparseLU solve failed";
  }
  return res;
}

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt) {
  SolveResult res;

  // Print header if verbose
  if (opt.verbose) {
    const char *method = opt.use_direct     ? "SparseLU"
                         : opt.use_bicgstab ? "BiCGSTAB"
                                            : "CG";
    const char *precond = (!opt.use_direct && opt.use_ilut) ? "+ILUT" : "";
    std::cerr << "Solver: Starting " << method << precond
              << " solve (N=" << A.rows() << ", nnz=" << A.nonZeros() << ")\n";
    std::cerr << std::fixed << std::setprecision(6);
  }

  if (opt.use_direct) {
    res = solve_direct(A, b, opt.verbose);
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
                  << " failed — falling back to SparseLU direct solver\n";
      }
      res = solve_direct(A, b, opt.verbose);
      res.method = failed_method + "->SparseLU";
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
                  << " failed — falling back to SparseLU direct solver\n";
      }
      res = solve_direct(A, b, opt.verbose);
      res.method = failed_method + "->SparseLU";
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

} // namespace edgefem
