#pragma once

#include <Eigen/Sparse>
#include <complex>
#include <unordered_set>

#include "vectorem/bc.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"

namespace vectorem {

using VecC = Eigen::VectorXcd;
using SpMatC = Eigen::SparseMatrix<std::complex<double>>;

struct MaxwellParams {
  double omega = 0.0;
  std::complex<double> eps_r = 1.0;
  std::complex<double> mu_r = 1.0;
  // Conductivity term for PML coordinate stretching (sigma)
  double pml_sigma = 0.0;
  // Physical region tags that should be treated as PML
  std::unordered_set<int> pml_regions;
  // Enable simple first-order absorbing boundary condition on outer faces
  bool use_abc = false;
};

struct MaxwellAssembly {
  SpMatC A;
  VecC b;
};

struct LumpedPort {
  int edge_id;
  double Z0 = 50.0;
};

MaxwellAssembly
assemble_maxwell(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                 const std::vector<LumpedPort> &ports,
                 int active_port_idx = -1); // -1 for no active port

Eigen::MatrixXcd
calculate_sparams(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                  const std::vector<LumpedPort> &ports);

} // namespace vectorem
