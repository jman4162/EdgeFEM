#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "vectorem/bc.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"
#include "vectorem/ports/wave_port.hpp"

namespace vectorem {

using VecC = Eigen::VectorXcd;
using SpMatC = Eigen::SparseMatrix<std::complex<double>>;

struct PMLRegionSpec {
  Eigen::Vector3d sigma_max = Eigen::Vector3d::Zero();
  Eigen::Vector3d thickness = Eigen::Vector3d::Zero();
  double grading_order = 3.0;
};

struct PMLDiagnostic {
  int region_tag = 0;
  Eigen::Vector3d sigma_max = Eigen::Vector3d::Zero();
  Eigen::Vector3d thickness = Eigen::Vector3d::Zero();
  Eigen::Vector3d reflection_est = Eigen::Vector3d::Ones();
};

struct MaxwellParams {
  double omega = 0.0;
  std::complex<double> eps_r = 1.0;
  std::complex<double> mu_r = 1.0;
  // Conductivity term for PML coordinate stretching (sigma)
  double pml_sigma = 0.0;
  // Physical region tags that should be treated as PML
  std::unordered_set<int> pml_regions;
  // Optional tensor PML specifications per physical region
  std::unordered_map<int, PMLRegionSpec> pml_tensor_regions;
  bool enforce_pml_heuristics = true;
  // Enable simple first-order absorbing boundary condition on outer faces
  bool use_abc = false;
};

struct MaxwellAssembly {
  SpMatC A;
  VecC b;
  std::vector<PMLDiagnostic> diagnostics;
};

MaxwellAssembly
assemble_maxwell(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                 const std::vector<WavePort> &ports,
                 int active_port_idx = -1); // -1 for no active port

Eigen::MatrixXcd
calculate_sparams(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                  const std::vector<WavePort> &ports);

} // namespace vectorem
