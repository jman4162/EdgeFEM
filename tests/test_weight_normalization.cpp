// Test S-parameters with different weight normalizations
// Key insight: V = w^H * x = A * ||w||² for pure mode with amplitude A
// For proper S-param extraction, we need V = sqrt(Z0) when A corresponds to unit power
// This requires ||w||² = sqrt(Z0) / A_unit, where A_unit = sqrt(2) (unit power normalization)
// So ||w||² = sqrt(Z0/2)

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

  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  double w1_norm_sq = wp1.weights.squaredNorm();
  double w2_norm_sq = wp2.weights.squaredNorm();
  double Z0 = std::real(wp1.mode.Z0);

  std::cout << "=== Port Weight Analysis ===" << std::endl;
  std::cout << "||w1||² = " << w1_norm_sq << std::endl;
  std::cout << "||w2||² = " << w2_norm_sq << std::endl;
  std::cout << "Z0 = " << Z0 << std::endl;
  std::cout << "sqrt(Z0) = " << std::sqrt(Z0) << std::endl;
  std::cout << "sqrt(Z0/2) = " << std::sqrt(Z0/2) << std::endl;
  std::cout << std::endl;

  MaxwellParams p;
  p.omega = omega;

  // Original (unnormalized)
  std::cout << "=== Original (||w||² = " << w1_norm_sq << ") ===" << std::endl;
  std::vector<WavePort> ports_orig{wp1, wp2};
  auto S_orig = calculate_sparams(mesh, p, bc, ports_orig);
  std::cout << "S11 = " << S_orig(0,0) << " (|S11| = " << std::abs(S_orig(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_orig(1,0) << " (|S21| = " << std::abs(S_orig(1,0)) << ")" << std::endl;

  // Normalize so ||w||² = sqrt(Z0)
  double target1 = std::sqrt(Z0);
  double scale1 = std::sqrt(target1 / w1_norm_sq);
  std::cout << "\n=== Normalized to ||w||² = sqrt(Z0) = " << target1 << " ===" << std::endl;
  std::cout << "Scale factor = " << scale1 << std::endl;

  WavePort wp1_norm1 = wp1, wp2_norm1 = wp2;
  wp1_norm1.weights *= scale1;
  wp2_norm1.weights *= scale1;

  std::vector<WavePort> ports_norm1{wp1_norm1, wp2_norm1};
  auto S_norm1 = calculate_sparams(mesh, p, bc, ports_norm1);
  std::cout << "S11 = " << S_norm1(0,0) << " (|S11| = " << std::abs(S_norm1(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_norm1(1,0) << " (|S21| = " << std::abs(S_norm1(1,0)) << ")" << std::endl;

  // Normalize so ||w||² = 1
  double target2 = 1.0;
  double scale2 = std::sqrt(target2 / w1_norm_sq);
  std::cout << "\n=== Normalized to ||w||² = 1 ===" << std::endl;
  std::cout << "Scale factor = " << scale2 << std::endl;

  WavePort wp1_norm2 = wp1, wp2_norm2 = wp2;
  wp1_norm2.weights *= scale2;
  wp2_norm2.weights *= scale2;

  std::vector<WavePort> ports_norm2{wp1_norm2, wp2_norm2};
  auto S_norm2 = calculate_sparams(mesh, p, bc, ports_norm2);
  std::cout << "S11 = " << S_norm2(0,0) << " (|S11| = " << std::abs(S_norm2(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_norm2(1,0) << " (|S21| = " << std::abs(S_norm2(1,0)) << ")" << std::endl;

  // Normalize so ||w||² = Z0
  double target3 = Z0;
  double scale3 = std::sqrt(target3 / w1_norm_sq);
  std::cout << "\n=== Normalized to ||w||² = Z0 = " << target3 << " ===" << std::endl;
  std::cout << "Scale factor = " << scale3 << std::endl;

  WavePort wp1_norm3 = wp1, wp2_norm3 = wp2;
  wp1_norm3.weights *= scale3;
  wp2_norm3.weights *= scale3;

  std::vector<WavePort> ports_norm3{wp1_norm3, wp2_norm3};
  auto S_norm3 = calculate_sparams(mesh, p, bc, ports_norm3);
  std::cout << "S11 = " << S_norm3(0,0) << " (|S11| = " << std::abs(S_norm3(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_norm3(1,0) << " (|S21| = " << std::abs(S_norm3(1,0)) << ")" << std::endl;

  // Try a range of normalizations to see the trend
  std::cout << "\n=== Sweep of ||w||² targets ===" << std::endl;
  std::cout << "Target        |S11|     |S21|     Passivity" << std::endl;

  for (double target : {0.1, 1.0, 10.0, 22.3, 100.0, 499.0, 1000.0, 5000.0, 38440.0}) {
    double s = std::sqrt(target / w1_norm_sq);
    WavePort wp1_t = wp1, wp2_t = wp2;
    wp1_t.weights *= s;
    wp2_t.weights *= s;

    std::vector<WavePort> ports_t{wp1_t, wp2_t};
    auto S_t = calculate_sparams(mesh, p, bc, ports_t);

    double pass = std::norm(S_t(0,0)) + std::norm(S_t(1,0));
    std::cout << std::setw(12) << target << "  "
              << std::setw(8) << std::abs(S_t(0,0)) << "  "
              << std::setw(8) << std::abs(S_t(1,0)) << "  "
              << (pass <= 1.01 ? "PASS" : "FAIL") << std::endl;
  }

  // Expected
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::cout << "\nExpected: |S11| = 0, |S21| = 1, phase(S21) = " << -beta*L*180/M_PI << "°" << std::endl;

  return 0;
}
