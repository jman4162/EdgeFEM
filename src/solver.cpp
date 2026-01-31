#include "vectorem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>

namespace vectorem {

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt) {
  SolveResult res;
  if (opt.use_bicgstab) {
    Eigen::BiCGSTAB<SpMatC> solver;
    solver.setTolerance(opt.tolerance);
    solver.setMaxIterations(opt.max_iterations);
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
      res.method = "BiCGSTAB";
      res.converged = false;
      res.error_message = "Matrix decomposition failed";
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
  return res;
}

} // namespace vectorem
