// Test S-parameters at different frequencies to understand FEM behavior
// Lower frequencies should work better if the issue is K >> k0²M

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  RectWaveguidePort dims{0.02286, 0.01016};
  double fc = c0 / (2 * dims.a);  // TE10 cutoff ~ 6.56 GHz

  std::cout << "TE10 cutoff frequency: " << fc/1e9 << " GHz" << std::endl;
  std::cout << "\nSweeping frequency from 7 GHz to 30 GHz:" << std::endl;
  std::cout << "freq(GHz)  k0²      K/(k0²M)  |S11|    |S21|    Phase(S21)  Status" << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  std::vector<double> freqs = {7e9, 8e9, 10e9, 15e9, 20e9, 25e9, 30e9};

  for (double freq : freqs) {
    if (freq <= fc) {
      std::cout << freq/1e9 << " GHz: Below cutoff, skipping" << std::endl;
      continue;
    }

    double omega = 2 * M_PI * freq;
    double k0 = omega / c0;
    double k0_sq = k0 * k0;

    // Get analytical mode parameters
    PortMode mode = solve_te10_mode(dims, freq);
    double beta = std::real(mode.beta);
    double Z0 = std::real(mode.Z0);

    // Build ports
    PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
    PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

    PortMode mode1 = solve_te10_mode(dims, freq);
    PortMode mode2 = solve_te10_mode(dims, freq);
    populate_te10_field(port1_surf, dims, mode1);
    populate_te10_field(port2_surf, dims, mode2);

    WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
    WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

    // Calculate K/(k0²M) ratio
    MaxwellParams p;
    p.omega = omega;
    auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);

    double K_avg = 0, M_est_avg = 0;
    int count = 0;
    for (int i = 0; i < (int)mesh.edges.size(); ++i) {
      if (bc.dirichlet_edges.count(i)) continue;
      K_avg += std::real(asmbl.A.coeff(i, i));
      count++;
    }
    K_avg /= count;

    // Estimate M from A = K - k0²M => M = (K - A)/k0² ≈ K/k0² (since A ≈ K for low freq)
    // Actually just use the formula: ratio = K_avg / (some estimated k0²M)
    // For simplicity, use K_avg directly as a proxy for the K contribution
    double km_ratio = K_avg / k0_sq;  // Rough estimate

    // Calculate S-parameters
    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    double s11_mag = std::abs(S(0,0));
    double s21_mag = std::abs(S(1,0));
    double s21_phase = std::arg(S(1,0)) * 180 / M_PI;

    // Expected S21 phase = -beta * L * 180/pi
    double L = 0.05;
    double s21_phase_expected = -beta * L * 180 / M_PI;
    while (s21_phase_expected < -180) s21_phase_expected += 360;

    // Check passivity
    double power_sum = s11_mag*s11_mag + s21_mag*s21_mag;
    std::string status = (power_sum <= 1.01) ? "PASS" : "FAIL";

    std::cout << std::setw(5) << freq/1e9 << "      "
              << std::setw(8) << k0_sq << "  "
              << std::setw(8) << km_ratio << "  "
              << std::setw(7) << s11_mag << "  "
              << std::setw(7) << s21_mag << "  "
              << std::setw(10) << s21_phase << "°  "
              << status << std::endl;
  }

  // Also test at very high frequency where mass term should dominate
  std::cout << "\n=== High Frequency Test (50-100 GHz) ===" << std::endl;
  std::vector<double> high_freqs = {50e9, 75e9, 100e9};

  for (double freq : high_freqs) {
    double omega = 2 * M_PI * freq;
    double k0 = omega / c0;
    double k0_sq = k0 * k0;

    PortMode mode = solve_te10_mode(dims, freq);

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

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    double s11_mag = std::abs(S(0,0));
    double s21_mag = std::abs(S(1,0));

    std::cout << freq/1e9 << " GHz: |S11|=" << s11_mag
              << ", |S21|=" << s21_mag
              << ", |S11|²+|S21|²=" << s11_mag*s11_mag + s21_mag*s21_mag
              << std::endl;
  }

  return 0;
}
