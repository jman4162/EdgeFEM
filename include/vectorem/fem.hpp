#pragma once

#include <Eigen/Sparse>
#include <complex>

#include "vectorem/bc.hpp"
#include "vectorem/mesh.hpp"

namespace vectorem {

using VecC = Eigen::VectorXcd;
using SpMatC = Eigen::SparseMatrix<std::complex<double>>;

struct ScalarHelmholtzParams {
  double k0 = 0.0; // wavenumber
  double eps_r = 1.0;
};

struct ScalarAssembly {
  SpMatC A;
  VecC b;
};

/// @warning STUB - Demo only. Returns identity matrix.
/// Use assemble_maxwell() for production simulations.
ScalarAssembly assemble_scalar_helmholtz(const Mesh &mesh,
                                         const ScalarHelmholtzParams &p,
                                         const BC &bc);

} // namespace vectorem
