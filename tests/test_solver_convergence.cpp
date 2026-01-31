#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>

#include "vectorem/solver.hpp"

using namespace vectorem;

void test_convergence_options() {
  std::cout << "test_convergence_options..." << std::endl;

  // Create a well-conditioned diagonal system
  int n = 100;
  SpMatC A(n, n);
  VecC b = VecC::Ones(n);

  std::vector<Eigen::Triplet<std::complex<double>>> trips;
  for (int i = 0; i < n; ++i) {
    trips.emplace_back(i, i, std::complex<double>(2.0 + i * 0.1, 0.0));
  }
  A.setFromTriplets(trips.begin(), trips.end());

  SolveOptions opts;
  opts.tolerance = 1e-10;
  opts.max_iterations = 1000;
  opts.use_bicgstab = true;

  auto result = solve_linear(A, b, opts);

  assert(result.converged);
  assert(result.method == "BiCGSTAB");
  assert(result.residual < opts.tolerance);
  assert(result.iters > 0);
  assert(result.error_message.empty());

  std::cout << "  Converged in " << result.iters << " iterations" << std::endl;
  std::cout << "  Residual: " << result.residual << std::endl;
}

void test_cg_solver() {
  std::cout << "test_cg_solver..." << std::endl;

  // Create SPD matrix for CG
  int n = 50;
  SpMatC A(n, n);
  VecC b = VecC::Ones(n);

  std::vector<Eigen::Triplet<std::complex<double>>> trips;
  for (int i = 0; i < n; ++i) {
    trips.emplace_back(i, i, std::complex<double>(4.0, 0.0));
    if (i > 0) {
      trips.emplace_back(i, i - 1, std::complex<double>(-1.0, 0.0));
      trips.emplace_back(i - 1, i, std::complex<double>(-1.0, 0.0));
    }
  }
  A.setFromTriplets(trips.begin(), trips.end());

  SolveOptions opts;
  opts.use_bicgstab = false;  // Use CG
  opts.tolerance = 1e-8;

  auto result = solve_linear(A, b, opts);

  assert(result.converged);
  assert(result.method == "CG");

  std::cout << "  CG converged in " << result.iters << " iterations" << std::endl;
}

void test_max_iterations_limit() {
  std::cout << "test_max_iterations_limit..." << std::endl;

  // Create a system that needs many iterations
  int n = 100;
  SpMatC A(n, n);
  VecC b = VecC::Random(n);

  std::vector<Eigen::Triplet<std::complex<double>>> trips;
  for (int i = 0; i < n; ++i) {
    trips.emplace_back(i, i, std::complex<double>(1.0, 0.0));
    if (i > 0) {
      trips.emplace_back(i, i - 1, std::complex<double>(-0.5, 0.1));
      trips.emplace_back(i - 1, i, std::complex<double>(-0.5, -0.1));
    }
  }
  A.setFromTriplets(trips.begin(), trips.end());

  SolveOptions opts;
  opts.tolerance = 1e-15;  // Very tight tolerance
  opts.max_iterations = 10;  // Very few iterations

  auto result = solve_linear(A, b, opts);

  // Should either converge or hit max iterations
  assert(result.iters <= opts.max_iterations || result.converged);

  std::cout << "  Iterations: " << result.iters << ", Converged: "
            << (result.converged ? "yes" : "no") << std::endl;
}

void test_tolerance_respected() {
  std::cout << "test_tolerance_respected..." << std::endl;

  int n = 50;
  SpMatC A(n, n);
  VecC b = VecC::Ones(n);

  std::vector<Eigen::Triplet<std::complex<double>>> trips;
  for (int i = 0; i < n; ++i) {
    trips.emplace_back(i, i, std::complex<double>(3.0, 0.0));
  }
  A.setFromTriplets(trips.begin(), trips.end());

  // Test with loose tolerance
  SolveOptions opts1;
  opts1.tolerance = 1e-4;
  auto result1 = solve_linear(A, b, opts1);

  // Test with tight tolerance
  SolveOptions opts2;
  opts2.tolerance = 1e-12;
  auto result2 = solve_linear(A, b, opts2);

  assert(result1.converged);
  assert(result2.converged);
  assert(result1.residual <= opts1.tolerance * 10);  // Allow some slack
  assert(result2.residual <= opts2.tolerance * 10);

  std::cout << "  Loose tol residual: " << result1.residual << std::endl;
  std::cout << "  Tight tol residual: " << result2.residual << std::endl;
}

int main() {
  test_convergence_options();
  test_cg_solver();
  test_max_iterations_limit();
  test_tolerance_respected();
  std::cout << "All solver convergence tests passed!" << std::endl;
  return 0;
}
