#pragma once

#include <functional>
#include <string>

#include "edgefem/fem.hpp"

namespace edgefem {

/// Callback function type for solver progress reporting.
/// @param iteration Current iteration number
/// @param residual Current residual norm
using SolverProgressCallback =
    std::function<void(int iteration, double residual)>;

struct SolveOptions {
  bool use_bicgstab = true;
  bool use_direct = false;  // Use SparseLU direct solver (more robust, uses more memory)
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

} // namespace edgefem
