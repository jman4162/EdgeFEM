#pragma once

#include <string>

#include "vectorem/fem.hpp"

namespace vectorem {

struct SolveOptions {
  bool use_bicgstab = true;
};

struct SolveResult {
  std::string method;
  int iters = 0;
  double residual = 0.0;
  VecC x;
};

SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt);

} // namespace vectorem
