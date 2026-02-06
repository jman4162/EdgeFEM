#include "edgefem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <iomanip>

namespace edgefem {

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt) {
  SolveResult res;

  // Print header if verbose
  if (opt.verbose) {
    std::cerr << "Solver: Starting " << (opt.use_bicgstab ? "BiCGSTAB" : "CG")
              << " solve (N=" << A.rows() << ", nnz=" << A.nonZeros()
              << ", tol=" << std::scientific << std::setprecision(2)
              << opt.tolerance << ", max_iter=" << opt.max_iterations << ")\n";
    std::cerr << std::fixed << std::setprecision(6);
  }

  if (opt.use_bicgstab) {
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
      return res;
    }
    res.x = solver.solve(b);
    res.method = "BiCGSTAB";
    res.iters = static_cast<int>(solver.iterations());
    res.residual = solver.error();
    res.converged = (solver.info() == Eigen::Success);
    if (!res.converged) {
      res.error_message = "Solver did not converge within max iterations";
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
      return res;
    }
    res.x = solver.solve(b);
    res.method = "CG";
    res.iters = static_cast<int>(solver.iterations());
    res.residual = solver.error();
    res.converged = (solver.info() == Eigen::Success);
    if (!res.converged) {
      res.error_message = "Solver did not converge within max iterations";
    }
  }

  // Print summary if verbose
  if (opt.verbose) {
    std::cerr << "Solver: " << res.method << " "
              << (res.converged ? "CONVERGED" : "FAILED")
              << " in " << res.iters << " iterations"
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
