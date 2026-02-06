// Test S-parameters with manual port sign correction
// The issue: both ports have +z normals, but port 1 should have -z
// This causes the weight signs to be wrong for S-parameter extraction
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

  std::cout << "=== Original (Both ports +z normal) ===" << std::endl;
  MaxwellParams p;
  p.omega = omega;
  std::vector<WavePort> ports_orig{wp1, wp2};
  auto S_orig = calculate_sparams(mesh, p, bc, ports_orig);
  std::cout << "S11 = " << S_orig(0,0) << " (|S11| = " << std::abs(S_orig(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_orig(1,0) << " (|S21| = " << std::abs(S_orig(1,0)) << ")" << std::endl;

  // Try flipping port 1 weights (should have -z normal)
  std::cout << "\n=== Flipped Port 1 (simulate -z normal) ===" << std::endl;
  WavePort wp1_flip = wp1;
  wp1_flip.weights *= -1.0;

  std::vector<WavePort> ports_flip1{wp1_flip, wp2};
  auto S_flip1 = calculate_sparams(mesh, p, bc, ports_flip1);
  std::cout << "S11 = " << S_flip1(0,0) << " (|S11| = " << std::abs(S_flip1(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_flip1(1,0) << " (|S21| = " << std::abs(S_flip1(1,0)) << ")" << std::endl;

  // Try flipping port 2 weights instead
  std::cout << "\n=== Flipped Port 2 ===" << std::endl;
  WavePort wp2_flip = wp2;
  wp2_flip.weights *= -1.0;

  std::vector<WavePort> ports_flip2{wp1, wp2_flip};
  auto S_flip2 = calculate_sparams(mesh, p, bc, ports_flip2);
  std::cout << "S11 = " << S_flip2(0,0) << " (|S11| = " << std::abs(S_flip2(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_flip2(1,0) << " (|S21| = " << std::abs(S_flip2(1,0)) << ")" << std::endl;

  // Try flipping both
  std::cout << "\n=== Flipped Both ===" << std::endl;
  std::vector<WavePort> ports_flip_both{wp1_flip, wp2_flip};
  auto S_flip_both = calculate_sparams(mesh, p, bc, ports_flip_both);
  std::cout << "S11 = " << S_flip_both(0,0) << " (|S11| = " << std::abs(S_flip_both(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_flip_both(1,0) << " (|S21| = " << std::abs(S_flip_both(1,0)) << ")" << std::endl;

  // Expected
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::complex<double> S21_expected = std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\n=== Expected (analytical) ===" << std::endl;
  std::cout << "S11 = 0" << std::endl;
  std::cout << "S21 = " << S21_expected << " (|S21| = 1, phase = " << std::arg(S21_expected)*180/M_PI << "°)" << std::endl;

  // Check passivity for each case
  std::cout << "\n=== Passivity Check ===" << std::endl;
  std::cout << "Original: |S11|² + |S21|² = " << std::norm(S_orig(0,0)) + std::norm(S_orig(1,0)) << std::endl;
  std::cout << "Flip1:    |S11|² + |S21|² = " << std::norm(S_flip1(0,0)) + std::norm(S_flip1(1,0)) << std::endl;
  std::cout << "Flip2:    |S11|² + |S21|² = " << std::norm(S_flip2(0,0)) + std::norm(S_flip2(1,0)) << std::endl;
  std::cout << "FlipBoth: |S11|² + |S21|² = " << std::norm(S_flip_both(0,0)) + std::norm(S_flip_both(1,0)) << std::endl;

  return 0;
}
