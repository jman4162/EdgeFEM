#pragma once

#include <Eigen/Sparse>
#include <complex>

#include "edgefem/bc.hpp"
#include "edgefem/mesh.hpp"

namespace edgefem {

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
/// Placeholder, not a Helmholtz assembly: returns the identity matrix with
/// Dirichlet rows set, and ignores p (k0, eps_r). Kept for the scalar demo
/// (EDGEFEM_BUILD_SCALAR, "scalar prototype"); do not use for physics.
ScalarAssembly assemble_scalar_helmholtz(const Mesh &mesh,
                                         const ScalarHelmholtzParams &p,
                                         const BC &bc);

} // namespace edgefem
