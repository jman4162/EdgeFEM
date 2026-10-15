#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "edgefem/bc.hpp"
#include "edgefem/materials/dispersive.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/periodic.hpp"
#include "edgefem/ports/port_eigensolve.hpp"
#include "edgefem/ports/wave_port.hpp"
#include "edgefem/solver.hpp"

namespace edgefem {

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
  None,           ///< No port ABC (default)
  Beta,           ///< j*beta (propagation constant)
  BetaNorm,       ///< j*beta/k0 (dimensionless form)
  ImpedanceMatch, ///< j*beta*sqrt(Z0/eta0) (impedance matched)
  ModalAdmittance ///< j*omega*eps0/Z0 (modal admittance based)
};

struct MaxwellParams {
  double omega = 0.0;
  std::complex<double> eps_r = 1.0;
  std::complex<double> mu_r = 1.0;

  // Region-based material assignment (takes precedence over global eps_r/mu_r)
  // Maps physical group tag to complex relative permittivity/permeability.
  // Elements with unassigned tags use the global eps_r/mu_r values above.
  std::unordered_map<int, std::complex<double>> eps_r_regions;
  std::unordered_map<int, std::complex<double>> mu_r_regions;

  // Dispersive material models per region (takes precedence over eps_r_regions)
  // When a region has a dispersive model, eps is evaluated at the current
  // omega. Maps physical group tag to dispersive material model.
  std::unordered_map<int, std::shared_ptr<materials::DispersiveMaterial>>
      eps_models;
  std::unordered_map<int, std::shared_ptr<materials::DispersiveMaterial>>
      mu_models;

  // Conductivity term for PML coordinate stretching (sigma)
  double pml_sigma = 0.0;
  // Physical region tags that should be treated as PML
  std::unordered_set<int> pml_regions;
  // Optional tensor PML specifications per physical region
  std::unordered_map<int, PMLRegionSpec> pml_tensor_regions;
  bool enforce_pml_heuristics = true;
  // Enable simple first-order absorbing boundary condition on outer faces
  bool use_abc = false;
  // Surface tags to apply ABC to. When non-empty, ABC is only applied to
  // triangles matching these physical tags. When empty, ABC applies to ALL
  // boundary triangles (backward-compatible default).
  std::unordered_set<int> abc_surface_tags;
  // Enable ABC damping at port edges (j*beta added to diagonal)
  // This improves wave absorption at ports for better S-parameter accuracy
  bool use_port_abc = false;
  // Port ABC formulation type (only used if use_port_abc = true)
  PortABCType port_abc_type = PortABCType::Beta;
  // Unused: no code path reads it (gap ports need no weight scaling). Kept
  // so existing callers that set it still compile.
  double port_weight_scale = 1.0;
  // Scaling of the port ABC operator jβ·M_s in calculate_sparams_eigenmode.
  // The correct value for the assembled surface mass matrix is 1.0; the
  // knob remains only as a diagnostic. (Before the assembled operator, a
  // diagonal-lumped approximation required an empirical 0.5 here.)
  double port_abc_scale = 1.0;
  // Use direct eigenmode excitation (more accurate, but slower)
  // When enabled, ports are excited using 3D FEM eigenvector instead of
  // analytical mode weights. Requires calling compute_te_eigenvector() first.
  bool use_eigenmode_excitation = false;

  /// Get effective permittivity for a given physical region tag
  /// (non-dispersive). Returns region-specific value if defined, otherwise
  /// global eps_r. Note: For dispersive materials, use get_eps_r(phys_tag,
  /// omega) instead.
  std::complex<double> get_eps_r(int phys_tag) const {
    auto it = eps_r_regions.find(phys_tag);
    return (it != eps_r_regions.end()) ? it->second : eps_r;
  }

  /// Get effective permittivity for a given physical region tag with frequency.
  /// Priority: dispersive model > static region value > global eps_r.
  /// @param phys_tag Physical group tag
  /// @param omega Angular frequency (rad/s) for dispersive evaluation
  std::complex<double> get_eps_r(int phys_tag, double omega) const {
    // First check for dispersive model
    auto model_it = eps_models.find(phys_tag);
    if (model_it != eps_models.end() && model_it->second) {
      return model_it->second->eval_eps(omega);
    }
    // Fall back to static region value or global
    auto it = eps_r_regions.find(phys_tag);
    return (it != eps_r_regions.end()) ? it->second : eps_r;
  }

  /// Get effective permeability for a given physical region tag
  /// (non-dispersive). Returns region-specific value if defined, otherwise
  /// global mu_r. Note: For dispersive materials, use get_mu_r(phys_tag, omega)
  /// instead.
  std::complex<double> get_mu_r(int phys_tag) const {
    auto it = mu_r_regions.find(phys_tag);
    return (it != mu_r_regions.end()) ? it->second : mu_r;
  }

  /// Get effective permeability for a given physical region tag with frequency.
  /// Priority: dispersive model > static region value > global mu_r.
  /// @param phys_tag Physical group tag
  /// @param omega Angular frequency (rad/s) for dispersive evaluation
  std::complex<double> get_mu_r(int phys_tag, double omega) const {
    // First check for dispersive model
    auto model_it = mu_models.find(phys_tag);
    if (model_it != mu_models.end() && model_it->second) {
      return model_it->second->eval_mu(omega);
    }
    // Fall back to static region value or global
    auto it = mu_r_regions.find(phys_tag);
    return (it != mu_r_regions.end()) ? it->second : mu_r;
  }
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
MaxwellAssembly assemble_maxwell_periodic(const Mesh &mesh,
                                          const MaxwellParams &p, const BC &bc,
                                          const PeriodicBC &pbc,
                                          const std::vector<WavePort> &ports,
                                          int active_port_idx = -1);

// ---------------------------------------------------------------------------
// Gap-port (lumped) network terms used by assemble_maxwell, calculate_sparams
// and frequency_sweep. The derivation is in src/assemble_maxwell.cpp. Each
// port's `weights` are its voltage functional: V = sum_k conj(weights_k) *
// x(edges_k) in volts (build_lumped_port sets weights = w / W, with w the
// surface-integral weights and W the strip width). Dirichlet (PEC) edges are
// skipped in all of these functions.
// ---------------------------------------------------------------------------

/// Frequency-independent port load L = sum_p (eta0 / Z0_p) v_p v_p^H over
/// ports with a nonzero Z0 and one weight per edge. The system matrix gets
/// j k0 L, the load of a resistor Z0_p across each port gap.
SpMatC port_load_matrix(const std::vector<WavePort> &ports, const BC &bc,
                        int num_edges);

/// Right-hand side for a unit incident power wave (a = 1) at `port`: a
/// Norton source of current 2 / sqrt(Z0) across the gap,
/// b = j k0 eta0 (2 / sqrt(Z0)) v.
VecC port_excitation(const WavePort &port, const BC &bc, int num_edges,
                     double k0);

/// Port voltage V = v^H x.
std::complex<double> port_voltage(const WavePort &port, const BC &bc,
                                  const VecC &x);

/// S-parameter column for an excitation of port `active` with a = 1:
/// S(j, active) = V_j / sqrt(Z0_j) - delta(j, active).
void fill_sparam_column(Eigen::MatrixXcd &S, int active,
                        const std::vector<WavePort> &ports, const BC &bc,
                        const VecC &x);

/// N-port S-matrix with gap (lumped) ports: the loaded system is assembled
/// and factored once (direct solver) and solved for one right-hand side per
/// port. Ports with modal weights also give a passive S here, but modal
/// ports are verified through calculate_sparams_eigenmode.
Eigen::MatrixXcd calculate_sparams(const Mesh &mesh, const MaxwellParams &p,
                                   const BC &bc,
                                   const std::vector<WavePort> &ports,
                                   const SolveOptions &opts = SolveOptions());

/// S-matrix and the field solutions it was computed from.
struct SParamSolution {
  Eigen::MatrixXcd S;
  std::vector<VecC> x; // x[k]: solution with port k excited (a_k = 1); empty
                       // if that solve failed
};

/// calculate_sparams() that also returns the field solutions, so a caller
/// that needs both the S-matrix and the fields factors the system once.
SParamSolution
calculate_sparams_with_fields(const Mesh &mesh, const MaxwellParams &p,
                              const BC &bc, const std::vector<WavePort> &ports,
                              const SolveOptions &opts = SolveOptions());

/// Deprecated: returns the ports unchanged and prints a one-time warning to
/// stderr. The gap-port formulation needs no weight normalization; the
/// former scaling (w^H A^-1 w = Z0) is what made a lossless port look
/// matched. Kept for one release so callers do not break.
void normalize_port_weights(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, std::vector<WavePort> &ports,
                            const SolveOptions &opts = SolveOptions());

/// Calculate S-parameters with periodic boundary conditions.
/// Useful for unit cell simulations under Floquet excitation.
Eigen::MatrixXcd calculate_sparams_periodic(const Mesh &mesh,
                                            const MaxwellParams &p,
                                            const BC &bc, const PeriodicBC &pbc,
                                            const std::vector<WavePort> &ports);

/// Calculate S-parameters using direct eigenmode excitation.
/// Uses 3D FEM eigenvectors as port weights, and a first-order modal ABC
/// with the ASSEMBLED port surface mass matrix: A += jβ·M_s per port,
/// RHS = 2jβ·M_s·e_inc, with e^H M_s e = 1 normalization and M_s-weighted
/// modal extraction. Still a single-mode first-order ABC: higher-order and
/// evanescent content at the port is not absorbed.
///
/// The ports must have been built with build_wave_port_from_eigenvector()
/// to contain the eigenvector-based weights and mode parameters (beta, kc,
/// etc).
///
/// @param mesh Volume mesh
/// @param p Maxwell parameters (omega required; leave port_abc_scale = 1.0)
/// @param bc PEC boundary conditions
/// @param ports Wave ports with eigenvector-based weights
/// @return N×N S-parameter matrix
Eigen::MatrixXcd
calculate_sparams_eigenmode(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports);

} // namespace edgefem
