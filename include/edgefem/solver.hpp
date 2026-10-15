#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "edgefem/fem.hpp"

namespace edgefem {

/// Callback function type for solver progress reporting.
/// @param iteration Current iteration number
/// @param residual Current residual norm
using SolverProgressCallback =
    std::function<void(int iteration, double residual)>;

/// Sparse LU used for direct solves (SolveOptions::use_direct and the
/// iterative solvers' fallback).
enum class DirectBackend {
  Default,  ///< UmfPack when EdgeFEM was built with SuiteSparse, else SparseLU
  SparseLU, ///< Eigen SparseLU
  UmfPack,  ///< SuiteSparse UMFPACK; std::invalid_argument without SuiteSparse
};

struct SolveOptions {
  bool use_bicgstab = true;
  bool use_direct =
      false; // Use a sparse direct LU (more robust, uses more memory)

  /// Which sparse LU a direct solve uses.
  DirectBackend direct_backend = DirectBackend::Default;
  double tolerance = 1e-10;
  int max_iterations = 10000;

  /// Use ILUT preconditioner for iterative solvers (BiCGSTAB/CG).
  /// Significantly improves convergence for ill-conditioned systems.
  bool use_ilut = true;

  /// ILUT fill factor: ratio of allowed fill-in to original nnz per row.
  /// Higher values improve preconditioner quality but use more memory.
  double ilut_fill_factor = 10.0;

  /// ILUT drop tolerance: entries smaller than this are dropped.
  /// Smaller values improve quality but increase memory/time.
  double ilut_drop_tolerance = 1e-4;

  /// Auto-fallback to SparseLU if iterative solver fails to converge.
  bool auto_fallback = true;

  /// Enable verbose output to stderr (prints iteration count and residual)
  bool verbose = false;

  /// Optional callback for progress monitoring (called every N iterations)
  SolverProgressCallback progress_callback = nullptr;

  /// How often to call progress_callback (every N iterations)
  int progress_interval = 100;
};

struct SolveResult {
  std::string method;
  int iters = 0;
  double residual = 0.0;
  bool converged = false;
  std::string error_message;
  VecC x;
};

/// Solve a complex sparse linear system Ax = b.
/// @param A System matrix (sparse, complex)
/// @param b Right-hand side vector
/// @param opt Solver options (tolerance, max iterations, verbose mode)
/// @return SolveResult containing solution vector and convergence info
SolveResult solve_linear(const SpMatC &A, const VecC &b,
                         const SolveOptions &opt = SolveOptions());

/// Solve A x_k = b_k for several right-hand sides. With opt.use_direct the
/// matrix is factored once and every b_k reuses the factorization; the
/// iterative path calls solve_linear for each b_k. A result whose
/// factorization failed has converged = false and an empty x.
std::vector<SolveResult>
solve_linear_multi(const SpMatC &A, const std::vector<VecC> &rhs,
                   const SolveOptions &opt = SolveOptions());

/// Warn on stderr when a solve failed: "WARNING: Solver did not converge
/// in <context> (method=..., iters=..., residual=..., error: ...)".
/// Returns res.converged.
bool check_solver_convergence(const SolveResult &res,
                              const std::string &context);

/// True when EdgeFEM was built with SuiteSparse UMFPACK
/// (-DEDGEFEM_USE_SUITESPARSE=ON; the default build is without it).
bool umfpack_available();

/// Number of sparse LU numeric factorizations (any backend) run in this
/// process. Every DirectLU::factorize call counts one. A diagnostic for
/// tests that check how many factorizations an operation costs; it is not
/// reset.
std::uint64_t direct_factorization_count();

/// Sparse LU factorization for matrices that share one sparsity pattern:
/// analyze the pattern once, factor each matrix, then solve any number of
/// right-hand sides with the last factorization. frequency_sweep uses one
/// instance for all frequencies.
///
/// Backends: UMFPACK (SuiteSparse, when built with it) and Eigen SparseLU
/// with COLAMD ordering.
class DirectLU {
public:
  /// @throws std::invalid_argument for DirectBackend::UmfPack in a build
  ///         without SuiteSparse.
  explicit DirectLU(DirectBackend backend = DirectBackend::Default);
  ~DirectLU();
  DirectLU(const DirectLU &) = delete;
  DirectLU &operator=(const DirectLU &) = delete;

  /// "UMFPACK" or "SparseLU".
  const std::string &method() const;

  /// Symbolic analysis of A's pattern. Later factorize() calls must pass
  /// matrices with this pattern; factorize() runs it if it has not run.
  void analyze_pattern(const SpMatC &A);

  /// Numeric factorization of A. Returns false if it failed (for example
  /// a singular matrix); solve() then fails until a factorization succeeds.
  bool factorize(const SpMatC &A);

  /// x = A^-1 b with the last successful factorization. Returns false (and
  /// leaves x empty) if there is none or the solve failed.
  bool solve(const VecC &b, VecC &x) const;

  /// factorize(A), then one solve per right-hand side. Each result carries
  /// method(), the relative residual, and converged = false with an empty
  /// x if the factorization or that solve failed.
  std::vector<SolveResult> factorize_and_solve(const SpMatC &A,
                                               const std::vector<VecC> &rhs);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace edgefem
