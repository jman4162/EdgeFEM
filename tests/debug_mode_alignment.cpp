// Compare port weights to FEM eigenvector for TE10 mode
// The hypothesis: port weights w don't align with FEM mode eigenvector

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

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
  PortMode mode1 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);

  // Get expected kc² for TE10
  double kc_te10 = M_PI / dims.a;
  double kc_te10_sq = kc_te10 * kc_te10;

  std::cout << "=== Mode Analysis ===" << std::endl;
  std::cout << "Expected TE10 kc² = " << kc_te10_sq << std::endl;
  std::cout << "k0² = " << k0*k0 << std::endl;
  std::cout << "Ratio k0²/kc² = " << k0*k0/kc_te10_sq << " (>1 means propagating)" << std::endl;

  // Assemble K and M matrices to find FEM eigenvector for TE10
  MaxwellParams p_low;
  p_low.omega = 1e-10;
  auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);

  double f_test = 1e9;
  double k0_test = 2*M_PI*f_test/c0;
  MaxwellParams p_test;
  p_test.omega = 2*M_PI*f_test;
  auto asmbl_test = assemble_maxwell(mesh, p_test, bc, {}, -1);

  int n = asmbl_low.A.rows();

  // Extract free DOFs
  std::vector<int> free_to_orig;
  std::vector<int> orig_to_free(n, -1);
  for (int i = 0; i < n; ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      orig_to_free[i] = free_to_orig.size();
      free_to_orig.push_back(i);
    }
  }
  int n_free = free_to_orig.size();

  std::cout << "Free DOFs: " << n_free << std::endl;

  // Build K and M for free DOFs
  Eigen::MatrixXd K_free(n_free, n_free);
  Eigen::MatrixXd M_free(n_free, n_free);

  for (int i = 0; i < n_free; ++i) {
    for (int j = 0; j < n_free; ++j) {
      int oi = free_to_orig[i];
      int oj = free_to_orig[j];
      double K_ij = std::real(asmbl_low.A.coeff(oi, oj));
      double A_test_ij = std::real(asmbl_test.A.coeff(oi, oj));
      double M_ij = (K_ij - A_test_ij) / (k0_test * k0_test);
      K_free(i, j) = K_ij;
      M_free(i, j) = M_ij;
    }
  }

  // Solve generalized eigenvalue problem
  std::cout << "\nSolving eigenvalue problem..." << std::endl;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);

  if (ges.info() != Eigen::Success) {
    std::cout << "Eigenvalue solver failed!" << std::endl;
    return 1;
  }

  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // Find eigenvalue closest to TE10 kc²
  int te10_idx = -1;
  double min_diff = 1e100;
  for (int i = 0; i < n_free; ++i) {
    if (eigenvalues(i) > 100) {  // Skip null space
      double diff = std::abs(eigenvalues(i) - kc_te10_sq);
      if (diff < min_diff) {
        min_diff = diff;
        te10_idx = i;
      }
    }
  }

  std::cout << "TE10 eigenvalue index: " << te10_idx << std::endl;
  std::cout << "TE10 eigenvalue: " << eigenvalues(te10_idx) << std::endl;
  std::cout << "Expected: " << kc_te10_sq << std::endl;
  std::cout << "Error: " << min_diff/kc_te10_sq*100 << "%" << std::endl;

  // Get the eigenvector
  Eigen::VectorXd v_te10 = eigenvectors.col(te10_idx);

  // Create full-size eigenvector (zero for PEC edges)
  Eigen::VectorXcd v_te10_full = Eigen::VectorXcd::Zero(n);
  for (int i = 0; i < n_free; ++i) {
    v_te10_full(free_to_orig[i]) = v_te10(i);
  }

  // Extract eigenvector values at port 1 edges
  std::cout << "\n=== Port Weight vs FEM Eigenvector ===" << std::endl;

  Eigen::VectorXcd w_port(wp1.edges.size());
  Eigen::VectorXcd v_port(wp1.edges.size());

  int free_port_edges = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    w_port(i) = wp1.weights(i);
    v_port(i) = v_te10_full(e);
    if (!bc.dirichlet_edges.count(e)) free_port_edges++;
  }

  std::cout << "Port edges: " << wp1.edges.size() << " (" << free_port_edges << " free)" << std::endl;

  // Normalize both for comparison
  double w_norm = w_port.norm();
  double v_norm = v_port.norm();

  std::cout << "||w_port|| = " << w_norm << std::endl;
  std::cout << "||v_port|| = " << v_norm << std::endl;

  if (w_norm > 1e-10 && v_norm > 1e-10) {
    Eigen::VectorXcd w_normalized = w_port / w_norm;
    Eigen::VectorXcd v_normalized = v_port / v_norm;

    // Compute correlation (inner product of normalized vectors)
    std::complex<double> correlation = w_normalized.dot(v_normalized);
    std::cout << "\nCorrelation <w,v> = " << correlation << std::endl;
    std::cout << "|<w,v>| = " << std::abs(correlation) << " (1.0 = perfect alignment)" << std::endl;

    // Also try with conjugate
    std::complex<double> correlation_conj = w_normalized.conjugate().dot(v_normalized);
    std::cout << "Correlation <w*,v> = " << correlation_conj << std::endl;
    std::cout << "|<w*,v>| = " << std::abs(correlation_conj) << std::endl;

    // Show first few pairs
    std::cout << "\nFirst 15 edge comparisons (free edges only):" << std::endl;
    std::cout << "Edge\t\tw (norm)\t\tv (norm)\t\tratio" << std::endl;
    int shown = 0;
    for (size_t i = 0; i < wp1.edges.size() && shown < 15; ++i) {
      int e = wp1.edges[i];
      if (!bc.dirichlet_edges.count(e)) {
        std::complex<double> w_n = w_normalized(i);
        std::complex<double> v_n = v_normalized(i);
        std::complex<double> ratio = (std::abs(v_n) > 1e-10) ? w_n / v_n : 0.0;
        std::cout << e << "\t\t" << w_n << "\t" << v_n << "\t" << ratio << std::endl;
        shown++;
      }
    }
  }

  // Now test: what if we use the FEM eigenvector as port weights?
  std::cout << "\n=== Testing FEM Eigenvector as Port Weights ===" << std::endl;

  // Create a port with FEM eigenvector as weights
  WavePort wp_fem = wp1;

  // Scale eigenvector so ||w||² = sqrt(Z0) for proper S-param extraction
  double Z0 = std::real(wp1.mode.Z0);
  double target_norm_sq = std::sqrt(Z0);

  // Extract only port edge values from eigenvector
  Eigen::VectorXcd v_weights(wp1.edges.size());
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    v_weights(i) = v_te10_full(wp1.edges[i]);
  }

  double v_port_norm_sq = v_weights.squaredNorm();
  double scale_factor = std::sqrt(target_norm_sq / v_port_norm_sq);
  v_weights *= scale_factor;

  std::cout << "Target ||w||² = sqrt(Z0) = " << target_norm_sq << std::endl;
  std::cout << "Scaled ||v_weights||² = " << v_weights.squaredNorm() << std::endl;

  // Also create port 2 with scaled eigenvector (need to compute for port 2)
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port2_surf, dims, mode2);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  // For port 2, extract eigenvector values and scale
  Eigen::VectorXcd v2_weights(wp2.edges.size());
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    v2_weights(i) = v_te10_full(wp2.edges[i]);
  }
  double v2_port_norm_sq = v2_weights.squaredNorm();
  double scale_factor2 = std::sqrt(target_norm_sq / v2_port_norm_sq);
  v2_weights *= scale_factor2;

  // Create modified ports
  WavePort wp1_fem = wp1;
  wp1_fem.weights = v_weights;

  WavePort wp2_fem = wp2;
  wp2_fem.weights = v2_weights;

  // Solve with FEM eigenvector weights
  MaxwellParams p;
  p.omega = omega;

  std::vector<WavePort> ports_fem{wp1_fem, wp2_fem};
  auto asmbl_fem = assemble_maxwell(mesh, p, bc, ports_fem, 0);
  auto res_fem = solve_linear(asmbl_fem.A, asmbl_fem.b, {});

  // Compute S-parameters
  std::complex<double> V1_fem = 0, V2_fem = 0;
  for (size_t i = 0; i < wp1_fem.edges.size(); ++i) {
    int e = wp1_fem.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1_fem += std::conj(wp1_fem.weights(i)) * res_fem.x(e);
    }
  }
  for (size_t i = 0; i < wp2_fem.edges.size(); ++i) {
    int e = wp2_fem.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2_fem += std::conj(wp2_fem.weights(i)) * res_fem.x(e);
    }
  }

  std::complex<double> V_inc_fem = std::sqrt(wp1.mode.Z0);
  std::complex<double> S11_fem = (V1_fem - V_inc_fem) / V_inc_fem;
  std::complex<double> S21_fem = V2_fem / V_inc_fem;

  std::cout << "\nS-parameters with FEM eigenvector weights:" << std::endl;
  std::cout << "V_inc = " << V_inc_fem << std::endl;
  std::cout << "V1 = " << V1_fem << std::endl;
  std::cout << "V2 = " << V2_fem << std::endl;
  std::cout << "S11 = " << S11_fem << " (|S11| = " << std::abs(S11_fem) << ")" << std::endl;
  std::cout << "S21 = " << S21_fem << " (|S21| = " << std::abs(S21_fem) << ")" << std::endl;

  // Expected
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::cout << "\nExpected: S11 = 0, |S21| = 1, phase(S21) = " << -beta*L*180/M_PI << "°" << std::endl;

  // Compare with original weights (properly normalized)
  std::cout << "\n=== Original Weights (normalized to ||w||² = sqrt(Z0)) ===" << std::endl;

  double orig_norm_sq = wp1.weights.squaredNorm();
  double orig_scale = std::sqrt(target_norm_sq / orig_norm_sq);

  WavePort wp1_norm = wp1;
  wp1_norm.weights *= orig_scale;
  WavePort wp2_norm = wp2;
  wp2_norm.weights *= orig_scale;

  std::vector<WavePort> ports_norm{wp1_norm, wp2_norm};
  auto asmbl_norm = assemble_maxwell(mesh, p, bc, ports_norm, 0);
  auto res_norm = solve_linear(asmbl_norm.A, asmbl_norm.b, {});

  std::complex<double> V1_norm = 0, V2_norm = 0;
  for (size_t i = 0; i < wp1_norm.edges.size(); ++i) {
    int e = wp1_norm.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1_norm += std::conj(wp1_norm.weights(i)) * res_norm.x(e);
    }
  }
  for (size_t i = 0; i < wp2_norm.edges.size(); ++i) {
    int e = wp2_norm.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2_norm += std::conj(wp2_norm.weights(i)) * res_norm.x(e);
    }
  }

  std::complex<double> S11_norm = (V1_norm - V_inc_fem) / V_inc_fem;
  std::complex<double> S21_norm = V2_norm / V_inc_fem;

  std::cout << "S11 = " << S11_norm << " (|S11| = " << std::abs(S11_norm) << ")" << std::endl;
  std::cout << "S21 = " << S21_norm << " (|S21| = " << std::abs(S21_norm) << ")" << std::endl;

  return 0;
}
