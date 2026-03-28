#pragma once

#include <complex>
#include <vector>

#include "edgefem/fem.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"
#include "edgefem/solver.hpp"

namespace edgefem {

enum class SweepPolicy { Discrete, Reuse, Balanced };

struct SweepLog {
  std::vector<double> solve_time;        // wall-clock seconds per solve
  std::vector<int> iterations;           // iteration count per solve
  std::vector<std::complex<double>> s21; // second port response
};

// Runs a sweep of a toy two-resonator system over the given frequencies.
// Returns timing and iteration logs as well as the S21 response at each point.
SweepLog sweep_two_resonator(const std::vector<double> &freqs,
                             SweepPolicy policy);

/// Pre-assembled stiffness (K) and mass (M) matrices for frequency sweeps.
/// A(ω) = K - k₀² M where k₀ = ω/c₀.
/// Only valid for non-dispersive materials without PML.
struct KMMatrices {
  SpMatC K; // Curl-curl stiffness matrix (frequency-independent)
  SpMatC M; // Mass matrix (frequency-independent)

  /// Combine K and M at a given frequency: A = K - k₀² M
  SpMatC combine(double omega) const;
};

/// Assemble K and M separately for frequency sweep reuse.
/// This avoids redundant element-level integration across frequencies.
/// Only valid when materials are non-dispersive and no PML is used.
///
/// @param mesh Volume mesh
/// @param p Maxwell params (eps_r, mu_r used; omega ignored for K/M assembly)
/// @param bc PEC boundary conditions
/// @return KMMatrices with pre-assembled K and M
KMMatrices assemble_maxwell_km(const Mesh &mesh, const MaxwellParams &p,
                               const BC &bc);

/// Result of a frequency sweep.
struct SweepResult {
  std::vector<double> frequencies;          // Hz
  std::vector<Eigen::MatrixXcd> S_matrices; // S-matrix at each frequency
};

/// Perform an S-parameter frequency sweep with matrix reuse.
///
/// Pre-assembles K and M once, then at each frequency:
/// 1. Combines A = K - k₀²M
/// 2. Adds port loading and ABC
/// 3. Solves and extracts S-parameters
///
/// For SparseLU direct solver, reuses symbolic factorization across
/// frequencies. For iterative solver, uses previous solution as initial guess.
///
/// @param mesh Volume mesh
/// @param p Maxwell params (omega is overwritten per frequency)
/// @param bc PEC boundary conditions
/// @param ports Wave/lumped port definitions
/// @param frequencies Vector of frequencies in Hz
/// @param opts Solver options
/// @return SweepResult with S-matrix at each frequency
SweepResult frequency_sweep(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports,
                            const std::vector<double> &frequencies,
                            const SolveOptions &opts = SolveOptions());

} // namespace edgefem
