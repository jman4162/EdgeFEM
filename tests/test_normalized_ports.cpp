// Test S-parameters with port weight normalization
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

  // Build ports
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  std::cout << "=== Before Normalization ===" << std::endl;
  double w1_norm_sq = 0;
  for (int i = 0; i < wp1.weights.size(); ++i)
    w1_norm_sq += std::norm(wp1.weights(i));
  std::cout << "||w1||² = " << w1_norm_sq << std::endl;
  std::cout << "Z0 = " << wp1.mode.Z0 << std::endl;

  // Calculate S-params without normalization
  MaxwellParams p;
  p.omega = omega;
  std::vector<WavePort> ports_unnorm{wp1, wp2};
  auto S_unnorm = calculate_sparams(mesh, p, bc, ports_unnorm);

  std::cout << "\nS-parameters (unnormalized):" << std::endl;
  std::cout << "S11 = " << S_unnorm(0,0) << " (|S11| = " << std::abs(S_unnorm(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_unnorm(1,0) << " (|S21| = " << std::abs(S_unnorm(1,0)) << ")" << std::endl;

  // Apply normalization
  std::cout << "\n=== Applying normalize_port_weights() ===" << std::endl;
  std::vector<WavePort> ports_norm{wp1, wp2};
  normalize_port_weights(mesh, p, bc, ports_norm);

  double w1_norm_sq_after = 0;
  for (int i = 0; i < ports_norm[0].weights.size(); ++i)
    w1_norm_sq_after += std::norm(ports_norm[0].weights(i));
  std::cout << "\n||w1||² after normalization = " << w1_norm_sq_after << std::endl;
  std::cout << "Scale factor = " << std::sqrt(w1_norm_sq_after / w1_norm_sq) << std::endl;

  // Calculate S-params with normalization
  auto S_norm = calculate_sparams(mesh, p, bc, ports_norm);

  std::cout << "\nS-parameters (normalized):" << std::endl;
  std::cout << "S11 = " << S_norm(0,0) << " (|S11| = " << std::abs(S_norm(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_norm(1,0) << " (|S21| = " << std::abs(S_norm(1,0)) << ")" << std::endl;
  std::cout << "|S11|² + |S21|² = " << std::norm(S_norm(0,0)) + std::norm(S_norm(1,0)) << std::endl;

  // Expected analytical
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::complex<double> S21_expected = std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\nExpected:" << std::endl;
  std::cout << "S11 = 0 (matched)" << std::endl;
  std::cout << "S21 = " << S21_expected << " (|S21| = 1)" << std::endl;

  // Try frequency sweep with normalization
  std::cout << "\n=== Frequency Sweep with Normalization ===" << std::endl;
  std::cout << "freq(GHz)  |S11|     |S21|     Passivity" << std::endl;

  for (double f : {7e9, 8e9, 10e9, 15e9, 20e9, 25e9, 30e9}) {
    double fc = c0 / (2 * dims.a);
    if (f <= fc) continue;

    double w = 2 * M_PI * f;

    PortMode m1 = solve_te10_mode(dims, f);
    PortMode m2 = solve_te10_mode(dims, f);
    populate_te10_field(port1_surf, dims, m1);
    populate_te10_field(port2_surf, dims, m2);

    WavePort wport1 = build_wave_port(mesh, port1_surf, m1);
    WavePort wport2 = build_wave_port(mesh, port2_surf, m2);

    MaxwellParams params;
    params.omega = w;

    std::vector<WavePort> ports{wport1, wport2};
    normalize_port_weights(mesh, params, bc, ports);

    auto S = calculate_sparams(mesh, params, bc, ports);

    double s11 = std::abs(S(0,0));
    double s21 = std::abs(S(1,0));
    double power = s11*s11 + s21*s21;

    std::cout << f/1e9 << " GHz:   " << s11 << "   " << s21 << "   "
              << (power <= 1.01 ? "PASS" : "FAIL") << std::endl;
  }

  return 0;
}
