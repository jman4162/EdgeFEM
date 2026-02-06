// Debug FEM matrix conditioning and eigenvalue structure
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;

  RectWaveguidePort dims{0.02286, 0.01016};
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);
  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  MaxwellParams p;
  p.omega = omega;

  // Assemble system with port loading
  std::vector<WavePort> ports{wp1, wp2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);

  int n = asmbl.A.rows();
  int n_free = n - bc.dirichlet_edges.size();

  std::cout << "Matrix size: " << n << " x " << n << std::endl;
  std::cout << "Free DOFs: " << n_free << std::endl;
  std::cout << "PEC DOFs: " << bc.dirichlet_edges.size() << std::endl;

  // Extract free-DOF submatrix for eigenvalue analysis
  // First, create mapping from original to free indices
  std::vector<int> free_to_orig;
  std::vector<int> orig_to_free(n, -1);
  for (int i = 0; i < n; ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      orig_to_free[i] = free_to_orig.size();
      free_to_orig.push_back(i);
    }
  }

  // Build dense free-DOF matrix
  Eigen::MatrixXcd A_free(n_free, n_free);
  Eigen::VectorXcd b_free(n_free);
  for (int i = 0; i < n_free; ++i) {
    int oi = free_to_orig[i];
    b_free(i) = asmbl.b(oi);
    for (int j = 0; j < n_free; ++j) {
      int oj = free_to_orig[j];
      A_free(i, j) = asmbl.A.coeff(oi, oj);
    }
  }

  std::cout << "\n=== Matrix Analysis ===" << std::endl;

  // Check symmetry
  double sym_error = (A_free - A_free.adjoint()).norm() / A_free.norm();
  std::cout << "Symmetry error: " << sym_error << std::endl;

  // Compute some eigenvalues of the system matrix
  // This is expensive but informative
  std::cout << "\nComputing eigenvalues (this may take a moment)..."
            << std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A_free);

  if (es.info() != Eigen::Success) {
    std::cout << "Eigenvalue computation failed!" << std::endl;
    return 1;
  }

  auto eigenvalues = es.eigenvalues();

  // Sort and display
  std::vector<double> sorted_eigs(eigenvalues.data(),
                                  eigenvalues.data() + n_free);
  std::sort(sorted_eigs.begin(), sorted_eigs.end());

  std::cout << "\nSmallest 10 eigenvalues:" << std::endl;
  for (int i = 0; i < std::min(10, n_free); ++i) {
    std::cout << "  λ_" << i << " = " << sorted_eigs[i] << std::endl;
  }

  std::cout << "\nLargest 10 eigenvalues:" << std::endl;
  for (int i = std::max(0, n_free - 10); i < n_free; ++i) {
    std::cout << "  λ_" << i << " = " << sorted_eigs[i] << std::endl;
  }

  double lambda_min = sorted_eigs[0];
  double lambda_max = sorted_eigs[n_free - 1];
  std::cout << "\nCondition number: " << lambda_max / lambda_min << std::endl;

  // Count positive, negative, and near-zero eigenvalues
  int n_pos = 0, n_neg = 0, n_zero = 0;
  for (double ev : sorted_eigs) {
    if (ev > 1e-6)
      n_pos++;
    else if (ev < -1e-6)
      n_neg++;
    else
      n_zero++;
  }
  std::cout << "\nEigenvalue signs: " << n_pos << " positive, " << n_neg
            << " negative, " << n_zero << " near-zero" << std::endl;

  // Expected: For wave propagation at 10 GHz, we need some negative eigenvalues
  // (modes below cutoff) and some positive (modes above cutoff)
  // If all positive, system is in quasi-static regime
  if (n_neg == 0) {
    std::cout << "\nWARNING: No negative eigenvalues found!" << std::endl;
    std::cout << "This means the system is positive definite, indicating"
              << std::endl;
    std::cout << "no propagating modes at this frequency." << std::endl;
    std::cout << "Expected negative eigenvalues for modes with kc < k0."
              << std::endl;
  }

  // What k0² eigenvalue would we need for TE10 mode?
  double kc_te10 = M_PI / dims.a;
  double kc_te10_sq = kc_te10 * kc_te10;
  double k0_sq = k0 * k0;
  std::cout << "\nTE10 mode kc² = " << kc_te10_sq << std::endl;
  std::cout << "Operating k0² = " << k0_sq << std::endl;
  std::cout << "For FEM matrix A = K - k0²M:" << std::endl;
  std::cout << "  Need eigenvalue of (K,M) near " << k0_sq << " for resonance"
            << std::endl;
  std::cout << "  Equivalently, need A eigenvalue near 0 for k² = k0²"
            << std::endl;

  // Find eigenvalue closest to zero
  double closest_to_zero = 1e100;
  for (double ev : sorted_eigs) {
    if (std::abs(ev) < std::abs(closest_to_zero))
      closest_to_zero = ev;
  }
  std::cout << "Closest eigenvalue to zero: " << closest_to_zero << std::endl;

  // Solve and check solution
  std::cout << "\n=== Solution Analysis ===" << std::endl;
  auto res = solve_linear(asmbl.A, asmbl.b, {});

  // Verify A*x = b
  VecC residual = asmbl.A * res.x - asmbl.b;
  double rel_residual = residual.norm() / asmbl.b.norm();
  std::cout << "Relative residual |Ax-b|/|b|: " << rel_residual << std::endl;

  // Check solution against null space of A
  // If solution has large component in null space direction, it's contaminated
  std::cout << "\nComparing solution to eigenvector structure..." << std::endl;
  Eigen::VectorXcd x_free(n_free);
  for (int i = 0; i < n_free; ++i) {
    x_free(i) = res.x(free_to_orig[i]);
  }

  // Project solution onto smallest eigenvectors
  auto eigenvectors = es.eigenvectors();
  std::cout << "Solution projection onto smallest eigenvectors:" << std::endl;
  for (int i = 0; i < std::min(5, n_free); ++i) {
    std::complex<double> proj = eigenvectors.col(i).adjoint() * x_free;
    std::cout << "  <v_" << i << ", x> = " << std::abs(proj) << std::endl;
  }

  return 0;
}
