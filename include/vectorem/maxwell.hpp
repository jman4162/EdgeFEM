#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "vectorem/bc.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/periodic.hpp"
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

/// ABC formulation variants for port boundary treatment
enum class PortABCType {
  None,         ///< No port ABC (default)
  Beta,         ///< j*beta (propagation constant)
  BetaNorm,     ///< j*beta/k0 (dimensionless form)
  ImpedanceMatch,///< j*beta*sqrt(Z0/eta0) (impedance matched)
  ModalAdmittance///< j*omega*eps0/Z0 (modal admittance based)
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
  // Enable ABC damping at port edges (j*beta added to diagonal)
  // This improves wave absorption at ports for better S-parameter accuracy
  bool use_port_abc = false;
  // Port ABC formulation type (only used if use_port_abc = true)
  PortABCType port_abc_type = PortABCType::Beta;
  // Port weight scaling factor (1.0 = use automatic normalization)
  // Set to custom value if manual tuning is needed
  double port_weight_scale = 1.0;
  // ABC coefficient scaling factor for port boundaries
  // Optimal value is ~0.5 due to diagonal approximation of surface mass matrix
  double port_abc_scale = 0.5;
  // Use direct eigenmode excitation (more accurate, but slower)
  // When enabled, ports are excited using 3D FEM eigenvector instead of
  // analytical mode weights. Requires calling compute_te_eigenvector() first.
  bool use_eigenmode_excitation = false;
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

/// Assemble Maxwell system with periodic boundary conditions.
/// Enforces phase-shifted periodic constraints via direct elimination:
/// slave DOFs are eliminated, their contributions accumulated to master DOFs.
///
/// @param mesh Volume mesh with periodic surfaces
/// @param p Maxwell parameters (frequency, materials, PML)
/// @param bc PEC boundary conditions
/// @param pbc Periodic boundary conditions with edge pairs and phase shift
/// @param ports Wave port definitions
/// @param active_port_idx Index of active port (-1 for no excitation)
/// @return Assembled system with periodic constraints applied
MaxwellAssembly
assemble_maxwell_periodic(const Mesh &mesh, const MaxwellParams &p,
                          const BC &bc, const PeriodicBC &pbc,
                          const std::vector<WavePort> &ports,
                          int active_port_idx = -1);

Eigen::MatrixXcd
calculate_sparams(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                  const std::vector<WavePort> &ports);

/// Normalize port weights for self-consistent port loading.
/// Computes w^H A^-1 w and scales weights so that w^H A^-1 w = Z0/2,
/// which is required for proper S-parameter extraction.
///
/// This should be called before calculate_sparams() for accurate results.
/// Without normalization, the port admittance may not properly match the
/// FEM matrix scale, leading to incorrect S-parameters.
///
/// @param mesh Volume mesh
/// @param p Maxwell parameters (frequency, materials)
/// @param bc PEC boundary conditions
/// @param ports Wave ports to normalize (modified in place)
void normalize_port_weights(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, std::vector<WavePort> &ports);

/// Calculate S-parameters with periodic boundary conditions.
/// Useful for unit cell simulations under Floquet excitation.
Eigen::MatrixXcd
calculate_sparams_periodic(const Mesh &mesh, const MaxwellParams &p,
                           const BC &bc, const PeriodicBC &pbc,
                           const std::vector<WavePort> &ports);

/// Calculate S-parameters using direct eigenmode excitation.
/// This method uses 3D FEM eigenvectors instead of analytical port weights,
/// providing better accuracy (~5% error) than the standard port formulation.
///
/// The ports must have been built with build_wave_port_from_eigenvector()
/// to contain the eigenvector-based weights and mode parameters (beta, kc, etc).
///
/// @param mesh Volume mesh
/// @param p Maxwell parameters (omega required, port_abc_scale recommended = 0.5)
/// @param bc PEC boundary conditions
/// @param ports Wave ports with eigenvector-based weights
/// @return N×N S-parameter matrix
Eigen::MatrixXcd
calculate_sparams_eigenmode(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc,
                            const std::vector<WavePort> &ports);

} // namespace vectorem
