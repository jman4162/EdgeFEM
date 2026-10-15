#pragma once

#include <complex>
#include <cstddef>
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

/// Field solutions to keep from a frequency sweep. By default none are kept.
struct SweepFieldRequest {
  /// Positions in the frequency list whose fields are kept.
  std::vector<std::size_t> indices;
  /// Also keep the fields at the frequency with the smallest |S(0,0)| (the
  /// first one on a tie; frequencies with a failed solve are skipped). The
  /// sweep holds only the best fields found so far, so this costs at most
  /// one extra set of field vectors.
  bool min_reflection = false;
};

/// Field solutions kept at one sweep frequency.
struct SweepField {
  std::size_t index = 0; // position in SweepResult::frequencies
  std::vector<VecC> x;   // x[k]: solution with port k excited (a_k = 1)
};

/// Result of a frequency sweep.
struct SweepResult {
  std::vector<double> frequencies;          // Hz
  std::vector<Eigen::MatrixXcd> S_matrices; // S-matrix at each frequency
  std::vector<SweepField> fields;           // kept fields, ascending index
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
/// @param keep Field solutions to return (none by default). Kept fields come
///        from the same factorization as the S-matrix: no extra solve.
/// @return SweepResult with S-matrix at each frequency and the kept fields
/// @throws std::invalid_argument if a kept index is out of range
SweepResult frequency_sweep(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports,
                            const std::vector<double> &frequencies,
                            const SolveOptions &opts = SolveOptions(),
                            const SweepFieldRequest &keep = SweepFieldRequest());

} // namespace edgefem
