// Final verification test for eigenmode-based S-parameter calculation
// Target criteria: |S11| < 0.1, |S21| > 0.95, passivity ≤ 1.01

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "vectorem/maxwell.hpp"
#include "vectorem/solver.hpp"

using namespace vectorem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;

int main() {
  std::cout << std::setprecision(6);

  // Load mesh and set up BCs
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  // Simulation parameters
  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double L = 0.05;  // 50mm waveguide length

  // WR-90 waveguide dimensions
  RectWaveguidePort dims{0.02286, 0.01016};
  double kc = M_PI / dims.a;
  double kc_sq = kc * kc;
  double beta = std::sqrt(k0*k0 - kc_sq);
  double Z0_te10 = omega * mu0 / beta;

  std::cout << "=== Eigenmode S-Parameter Verification Test ===" << std::endl;
  std::cout << "WR-90 Waveguide, 50mm length, 10 GHz" << std::endl;
  std::cout << "TE10 beta = " << beta << " rad/m" << std::endl;
  std::cout << "TE10 Z0 = " << Z0_te10 << " ohms" << std::endl;
  std::cout << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Compute 3D FEM eigenvector for TE10 mode
  std::cout << "Computing 3D FEM eigenvector for TE10 mode..." << std::endl;
  Eigen::VectorXd v_te10 = compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);
  std::cout << "Eigenvector norm: " << v_te10.norm() << std::endl;

  // Create analytical mode parameters
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);

  // Build ports using eigenvector-based weights
  std::cout << "Building ports with eigenvector-based weights..." << std::endl;
  WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10, mode1, bc.dirichlet_edges);
  WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10, mode2, bc.dirichlet_edges);

  std::cout << "Port 1: " << wp1.edges.size() << " edges" << std::endl;
  std::cout << "Port 2: " << wp2.edges.size() << " edges" << std::endl;

  // Set up Maxwell parameters with optimal ABC scale
  MaxwellParams p;
  p.omega = omega;
  p.port_abc_scale = 0.5;  // Optimal value from tuning

  std::vector<WavePort> ports{wp1, wp2};

  // Calculate S-parameters using eigenmode method
  std::cout << "\nCalculating S-parameters..." << std::endl;
  auto S = calculate_sparams_eigenmode(mesh, p, bc, ports);

  // Results
  std::cout << "\n=== S-Parameter Results ===" << std::endl;
  std::cout << "S11 = " << S(0,0) << " (|S11| = " << std::abs(S(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S(1,0) << " (|S21| = " << std::abs(S(1,0))
            << ", phase = " << std::arg(S(1,0))*180/M_PI << "°)" << std::endl;
  std::cout << "S12 = " << S(0,1) << " (|S12| = " << std::abs(S(0,1)) << ")" << std::endl;
  std::cout << "S22 = " << S(1,1) << " (|S22| = " << std::abs(S(1,1)) << ")" << std::endl;

  // Expected values
  double expected_phase = -beta * L * 180 / M_PI;
  // Normalize to -180 to 180 range
  while (expected_phase < -180) expected_phase += 360;
  while (expected_phase > 180) expected_phase -= 360;

  std::cout << "\n=== Expected Values ===" << std::endl;
  std::cout << "S11 = 0" << std::endl;
  std::cout << "S21: |S21| = 1, phase = " << expected_phase << "°" << std::endl;

  // Verification
  std::cout << "\n=== Verification ===" << std::endl;

  double s11_mag = std::abs(S(0,0));
  double s21_mag = std::abs(S(1,0));
  double passivity = std::norm(S(0,0)) + std::norm(S(1,0));
  double phase_error = std::abs(std::arg(S(1,0))*180/M_PI - expected_phase);
  if (phase_error > 180) phase_error = 360 - phase_error;

  std::cout << "|S11| = " << s11_mag << " (target: < 0.1) "
            << (s11_mag < 0.1 ? "✓ PASS" : "✗ FAIL") << std::endl;
  std::cout << "|S21| = " << s21_mag << " (target: > 0.95) "
            << (s21_mag > 0.95 ? "✓ PASS" : "✗ FAIL") << std::endl;
  std::cout << "Passivity = " << passivity << " (target: ≤ 1.01) "
            << (passivity <= 1.01 ? "✓ PASS" : "✗ FAIL") << std::endl;
  std::cout << "Phase error = " << phase_error << "° (target: < 10°) "
            << (phase_error < 10 ? "✓ PASS" : "✗ FAIL") << std::endl;

  // Reciprocity check
  double s12_mag = std::abs(S(0,1));
  double s21_s12_diff = std::abs(s21_mag - s12_mag);
  std::cout << "|S12 - S21| = " << s21_s12_diff << " (target: < 0.05) "
            << (s21_s12_diff < 0.05 ? "✓ PASS" : "✗ FAIL") << std::endl;

  // Overall result
  bool all_pass = (s11_mag < 0.1) && (s21_mag > 0.95) && (passivity <= 1.01) && (phase_error < 10);

  std::cout << "\n=== OVERALL: " << (all_pass ? "PASS" : "FAIL") << " ===" << std::endl;

  return all_pass ? 0 : 1;
}
