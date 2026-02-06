// Test combining wave port with ABC to help absorption
// The hypothesis: wave port loading alone is too weak (0.8% of matrix norm)
// ABC adds imaginary damping to boundary edges that may help

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

  std::cout << "=== Baseline: No ABC ===" << std::endl;
  MaxwellParams p_noabc;
  p_noabc.omega = omega;
  p_noabc.use_abc = false;

  std::vector<WavePort> ports{wp1, wp2};
  auto S_noabc = calculate_sparams(mesh, p_noabc, bc, ports);
  std::cout << "S11 = " << S_noabc(0,0) << " (|S11| = " << std::abs(S_noabc(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_noabc(1,0) << " (|S21| = " << std::abs(S_noabc(1,0)) << ")" << std::endl;

  std::cout << "\n=== With ABC on boundary ===" << std::endl;
  MaxwellParams p_abc;
  p_abc.omega = omega;
  p_abc.use_abc = true;

  auto S_abc = calculate_sparams(mesh, p_abc, bc, ports);
  std::cout << "S11 = " << S_abc(0,0) << " (|S11| = " << std::abs(S_abc(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_abc(1,0) << " (|S21| = " << std::abs(S_abc(1,0)) << ")" << std::endl;

  // Try scaling up the port weights to make port loading stronger
  std::cout << "\n=== Increased Port Loading (10x weights) ===" << std::endl;
  WavePort wp1_10x = wp1, wp2_10x = wp2;
  wp1_10x.weights *= 10.0;
  wp2_10x.weights *= 10.0;

  std::vector<WavePort> ports_10x{wp1_10x, wp2_10x};
  auto S_10x = calculate_sparams(mesh, p_noabc, bc, ports_10x);
  std::cout << "S11 = " << S_10x(0,0) << " (|S11| = " << std::abs(S_10x(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_10x(1,0) << " (|S21| = " << std::abs(S_10x(1,0)) << ")" << std::endl;

  // Try even more scaling
  std::cout << "\n=== Increased Port Loading (100x weights) ===" << std::endl;
  WavePort wp1_100x = wp1, wp2_100x = wp2;
  wp1_100x.weights *= 100.0;
  wp2_100x.weights *= 100.0;

  std::vector<WavePort> ports_100x{wp1_100x, wp2_100x};
  auto S_100x = calculate_sparams(mesh, p_noabc, bc, ports_100x);
  std::cout << "S11 = " << S_100x(0,0) << " (|S11| = " << std::abs(S_100x(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_100x(1,0) << " (|S21| = " << std::abs(S_100x(1,0)) << ")" << std::endl;

  // Try a different approach: normalize weights so ||w||²/Z0 matches matrix diagonal
  // Matrix diagonal ~ 1000, so we want ||w||²/Z0 ~ 1000, i.e., ||w||² ~ 500000
  std::cout << "\n=== Matched Port Loading (||w||²/Z0 ~ A_ii) ===" << std::endl;
  double Z0 = std::real(wp1.mode.Z0);
  double A_diag = 1000.0;  // Typical diagonal value
  double target_w_sq = A_diag * Z0;  // ~ 500000
  double current_w_sq = wp1.weights.squaredNorm();  // ~ 38440
  double scale_match = std::sqrt(target_w_sq / current_w_sq);

  std::cout << "Target ||w||² = " << target_w_sq << std::endl;
  std::cout << "Scale factor = " << scale_match << std::endl;

  WavePort wp1_match = wp1, wp2_match = wp2;
  wp1_match.weights *= scale_match;
  wp2_match.weights *= scale_match;

  std::vector<WavePort> ports_match{wp1_match, wp2_match};
  auto S_match = calculate_sparams(mesh, p_noabc, bc, ports_match);
  std::cout << "S11 = " << S_match(0,0) << " (|S11| = " << std::abs(S_match(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_match(1,0) << " (|S21| = " << std::abs(S_match(1,0)) << ")" << std::endl;

  // Expected values
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::complex<double> S21_expected = std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\nExpected: S11 = 0, S21 = " << S21_expected << std::endl;
  std::cout << "Phase(S21_expected) = " << std::arg(S21_expected)*180/M_PI << "°" << std::endl;

  // Check passivity for different scalings
  std::cout << "\n=== Passivity Check ===" << std::endl;
  std::cout << "No ABC:  |S11|² + |S21|² = " << std::norm(S_noabc(0,0)) + std::norm(S_noabc(1,0)) << std::endl;
  std::cout << "ABC:     |S11|² + |S21|² = " << std::norm(S_abc(0,0)) + std::norm(S_abc(1,0)) << std::endl;
  std::cout << "10x:     |S11|² + |S21|² = " << std::norm(S_10x(0,0)) + std::norm(S_10x(1,0)) << std::endl;
  std::cout << "100x:    |S11|² + |S21|² = " << std::norm(S_100x(0,0)) + std::norm(S_100x(1,0)) << std::endl;
  std::cout << "Matched: |S11|² + |S21|² = " << std::norm(S_match(0,0)) + std::norm(S_match(1,0)) << std::endl;

  return 0;
}
