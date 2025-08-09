#pragma once

#include <Eigen/Sparse>
#include <complex>

#include "vectorem/bc.hpp"
#include "vectorem/mesh.hpp"

namespace vectorem {

using VecC = Eigen::VectorXcd;
using SpMatC = Eigen::SparseMatrix<std::complex<double>>;

struct MaxwellParams {
  double omega = 0.0;
  double eps_r = 1.0;
  double mu_r = 1.0;
};

struct MaxwellAssembly {
  SpMatC A;
  VecC b;
};

MaxwellAssembly assemble_maxwell(const Mesh &mesh,
                                 const MaxwellParams &p,
                                 const BC &bc);

} // namespace vectorem
