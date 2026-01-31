#pragma once

#include <string>

#include "vectorem/fem.hpp"

namespace vectorem {

struct SolveOptions {
  bool use_bicgstab = true;
  double tolerance = 1e-10;
  int max_iterations = 10000;
};

struct SolveResult {
  std::string method;
  int iters = 0;
  double residual = 0.0;
  bool converged = false;
  std::string error_message;
  VecC x;
};

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt = SolveOptions());

} // namespace vectorem
