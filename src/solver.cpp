#include "vectorem/solver.hpp"

#include <Eigen/IterativeLinearSolvers>

namespace vectorem {

SolveResult solve_linear(const SpMatC &A, const VecC &b, const SolveOptions &opt) {
    SolveResult res;
    if (opt.use_bicgstab) {
        Eigen::BiCGSTAB<SpMatC> solver;
        solver.compute(A);
        res.x = solver.solve(b);
        res.method = "BiCGSTAB";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
    } else {
        Eigen::ConjugateGradient<SpMatC> solver;
        solver.compute(A);
        res.x = solver.solve(b);
        res.method = "CG";
        res.iters = static_cast<int>(solver.iterations());
        res.residual = solver.error();
    }
    return res;
}

} // namespace vectorem

