// Test the new eigenvector-based wave port formulation
// This should give much better S-parameter accuracy than analytical weights

#include <iostream>
#include <iomanip>
#include <cmath>
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

  RectWaveguidePort dims{0.02286, 0.01016};

  // Compute TE10 cutoff
  double kc_te10_sq = std::pow(M_PI / dims.a, 2);

  std::cout << "=== Eigenvector-Based Wave Port Test ===" << std::endl;
  std::cout << "Frequency: " << freq/1e9 << " GHz" << std::endl;
  std::cout << "Target kc² (TE10): " << kc_te10_sq << std::endl;

  // Compute 3D FEM eigenvector for TE10 mode
  std::cout << "\nComputing 3D FEM eigenvector..." << std::endl;
  Eigen::VectorXd v_te10 = compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_te10_sq);

  double v_norm = v_te10.norm();
  std::cout << "Eigenvector norm: " << v_norm << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Create analytical mode parameters
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);

  // Build ports with eigenvector-based weights
  WavePort wp1_eig = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10, mode1, bc.dirichlet_edges);
  WavePort wp2_eig = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10, mode2, bc.dirichlet_edges);

  std::cout << "\nPort 1 eigenvector weights ||w||² = " << wp1_eig.weights.squaredNorm() << std::endl;
  std::cout << "Port 2 eigenvector weights ||w||² = " << wp2_eig.weights.squaredNorm() << std::endl;
  std::cout << "Target ||w||² = sqrt(Z0) = " << std::sqrt(std::real(mode1.Z0)) << std::endl;

  // Also build analytical ports for comparison
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);
  WavePort wp1_ana = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2_ana = build_wave_port(mesh, port2_surf, mode2);

  // Normalize analytical weights to same scale
  double ana_norm_sq = wp1_ana.weights.squaredNorm();
  double target = std::sqrt(std::real(mode1.Z0));
  double ana_scale = std::sqrt(target / ana_norm_sq);
  wp1_ana.weights *= ana_scale;
  wp2_ana.weights *= ana_scale;

  std::cout << "Analytical weights (normalized) ||w||² = " << wp1_ana.weights.squaredNorm() << std::endl;

  // Test 1: Analytical weights without port ABC
  std::cout << "\n=== Test 1: Analytical Weights (no port ABC) ===" << std::endl;
  MaxwellParams p1;
  p1.omega = omega;
  p1.use_port_abc = false;

  std::vector<WavePort> ports_ana{wp1_ana, wp2_ana};
  auto S_ana = calculate_sparams(mesh, p1, bc, ports_ana);
  std::cout << "S11 = " << S_ana(0,0) << " (|S11| = " << std::abs(S_ana(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_ana(1,0) << " (|S21| = " << std::abs(S_ana(1,0)) << ")" << std::endl;

  // Test 2: Analytical weights with port ABC
  std::cout << "\n=== Test 2: Analytical Weights (with port ABC) ===" << std::endl;
  MaxwellParams p2;
  p2.omega = omega;
  p2.use_port_abc = true;

  auto S_ana_abc = calculate_sparams(mesh, p2, bc, ports_ana);
  std::cout << "S11 = " << S_ana_abc(0,0) << " (|S11| = " << std::abs(S_ana_abc(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_ana_abc(1,0) << " (|S21| = " << std::abs(S_ana_abc(1,0)) << ")" << std::endl;

  // Test 3: Eigenvector weights without port ABC
  std::cout << "\n=== Test 3: Eigenvector Weights (no port ABC) ===" << std::endl;
  std::vector<WavePort> ports_eig{wp1_eig, wp2_eig};
  auto S_eig = calculate_sparams(mesh, p1, bc, ports_eig);
  std::cout << "S11 = " << S_eig(0,0) << " (|S11| = " << std::abs(S_eig(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_eig(1,0) << " (|S21| = " << std::abs(S_eig(1,0)) << ")" << std::endl;

  // Test 4: Eigenvector weights with port ABC
  std::cout << "\n=== Test 4: Eigenvector Weights (with port ABC) ===" << std::endl;
  auto S_eig_abc = calculate_sparams(mesh, p2, bc, ports_eig);
  std::cout << "S11 = " << S_eig_abc(0,0) << " (|S11| = " << std::abs(S_eig_abc(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_eig_abc(1,0) << " (|S21| = " << std::abs(S_eig_abc(1,0)) << ")" << std::endl;

  // Expected values
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::complex<double> S21_expected = std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\n=== Expected ===" << std::endl;
  std::cout << "S11 = 0" << std::endl;
  std::cout << "S21 = " << S21_expected << " (|S21| = 1, phase = " << std::arg(S21_expected)*180/M_PI << "°)" << std::endl;

  // Passivity check
  std::cout << "\n=== Passivity Check (|S11|² + |S21|² ≤ 1) ===" << std::endl;
  std::cout << "Analytical:     " << std::norm(S_ana(0,0)) + std::norm(S_ana(1,0)) << std::endl;
  std::cout << "Ana + ABC:      " << std::norm(S_ana_abc(0,0)) + std::norm(S_ana_abc(1,0)) << std::endl;
  std::cout << "Eigenvector:    " << std::norm(S_eig(0,0)) + std::norm(S_eig(1,0)) << std::endl;
  std::cout << "Eig + ABC:      " << std::norm(S_eig_abc(0,0)) + std::norm(S_eig_abc(1,0)) << std::endl;

  // Check correlation between analytical and eigenvector weights
  std::cout << "\n=== Weight Correlation ===" << std::endl;
  std::complex<double> corr = wp1_ana.weights.dot(wp1_eig.weights);
  corr /= (wp1_ana.weights.norm() * wp1_eig.weights.norm());
  std::cout << "|<w_ana, w_eig>| = " << std::abs(corr) << std::endl;

  return 0;
}
