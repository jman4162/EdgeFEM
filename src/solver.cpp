#include "edgefem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <iomanip>
#include <iostream>

namespace edgefem {

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt) {
  SolveResult res;

  // Print header if verbose
  if (opt.verbose) {
    const char *method = opt.use_direct ? "SparseLU"
                         : opt.use_bicgstab ? "BiCGSTAB" : "CG";
    std::cerr << "Solver: Starting " << method
              << " solve (N=" << A.rows() << ", nnz=" << A.nonZeros()
              << ")\n";
    std::cerr << std::fixed << std::setprecision(6);
  }

  if (opt.use_direct) {
    Eigen::SparseLU<SpMatC> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
      res.method = "SparseLU";
      res.converged = false;
      res.error_message = "SparseLU factorization failed";
      if (opt.verbose) {
        std::cerr << "Solver: SparseLU factorization FAILED\n";
      }
      return res;
    }
    res.x = solver.solve(b);
    res.method = "SparseLU";
    res.iters = 1;
    res.residual = (A * res.x - b).norm() / b.norm();
    res.converged = (solver.info() == Eigen::Success);
    if (!res.converged) {
      res.error_message = "SparseLU solve failed";
    }
  } else if (opt.use_bicgstab) {
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
