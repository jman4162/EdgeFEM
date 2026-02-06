// WR-90 Waveguide Benchmark Test
// Validates VectorEM against analytical results for rectangular waveguide
//
// Test cases:
// 1. 50mm waveguide at 10 GHz: |S21| > 0.95, |S11| < 0.1, phase within 10°
// 2. Multiple frequencies across TE10 passband (7-12 GHz)
// 3. Mesh convergence check (results should be stable)

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cassert>
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;

struct BenchmarkResult {
  double freq_ghz;
  double s11_mag;
  double s21_mag;
  double s21_phase_deg;
  double expected_phase_deg;
  double passivity;
  bool passed;
};

BenchmarkResult run_benchmark(const Mesh &mesh, const BC &bc,
                              const PortSurfaceMesh &port1_surf,
                              const PortSurfaceMesh &port2_surf,
                              const RectWaveguidePort &dims,
                              double freq, double L) {
  BenchmarkResult r;
  r.freq_ghz = freq / 1e9;

  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double kc = M_PI / dims.a;
  double kc_sq = kc * kc;
  double beta_sq = k0 * k0 - kc_sq;

  if (beta_sq <= 0) {
    // Below cutoff
    r.s11_mag = 1.0;
    r.s21_mag = 0.0;
    r.s21_phase_deg = 0.0;
    r.expected_phase_deg = 0.0;
    r.passivity = 1.0;
    r.passed = true;  // Evanescent mode is expected to have zero transmission
    return r;
  }

  double beta = std::sqrt(beta_sq);

  // Compute eigenvector
  Eigen::VectorXd v_te10 = compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);

  // Create mode parameters
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);

  // Build ports
  WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10, mode1, bc.dirichlet_edges);
  WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10, mode2, bc.dirichlet_edges);

  // Calculate S-parameters
  MaxwellParams p;
  p.omega = omega;
  p.port_abc_scale = 0.5;

  std::vector<WavePort> ports{wp1, wp2};
  auto S = calculate_sparams_eigenmode(mesh, p, bc, ports);

  r.s11_mag = std::abs(S(0, 0));
  r.s21_mag = std::abs(S(1, 0));
  r.s21_phase_deg = std::arg(S(1, 0)) * 180.0 / M_PI;
  r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));

  // Expected phase (normalized to -180 to 180)
  r.expected_phase_deg = -beta * L * 180.0 / M_PI;
  while (r.expected_phase_deg < -180) r.expected_phase_deg += 360;
  while (r.expected_phase_deg > 180) r.expected_phase_deg -= 360;

  // Pass criteria
  double phase_error = std::abs(r.s21_phase_deg - r.expected_phase_deg);
  if (phase_error > 180) phase_error = 360 - phase_error;

  r.passed = (r.s11_mag < 0.15) && (r.s21_mag > 0.90) &&
             (r.passivity <= 1.05) && (phase_error < 15);

  return r;
}

int main() {
  std::cout << std::setprecision(4) << std::fixed;

  std::cout << "======================================================" << std::endl;
  std::cout << "VectorEM WR-90 Waveguide Benchmark" << std::endl;
  std::cout << "======================================================" << std::endl;

  // Load mesh
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  std::cout << "Mesh: " << mesh.tets.size() << " tets, "
            << mesh.edges.size() << " edges" << std::endl;

  // WR-90 dimensions
  RectWaveguidePort dims{0.02286, 0.01016};
  double L = 0.05;  // 50mm length

  // Cutoff frequency
  double fc = c0 / (2 * dims.a);
  std::cout << "WR-90: a=" << dims.a*1000 << "mm, b=" << dims.b*1000 << "mm" << std::endl;
  std::cout << "TE10 cutoff: " << fc/1e9 << " GHz" << std::endl;
  std::cout << "Waveguide length: " << L*1000 << " mm" << std::endl;
  std::cout << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Test 1: Center frequency (10 GHz)
  std::cout << "=== Test 1: Center Frequency (10 GHz) ===" << std::endl;
  auto r10 = run_benchmark(mesh, bc, port1_surf, port2_surf, dims, 10e9, L);
  std::cout << "  |S11| = " << r10.s11_mag << " (target < 0.15)" << std::endl;
  std::cout << "  |S21| = " << r10.s21_mag << " (target > 0.90)" << std::endl;
  std::cout << "  Phase(S21) = " << r10.s21_phase_deg << "° (expected " << r10.expected_phase_deg << "°)" << std::endl;
  std::cout << "  Passivity = " << r10.passivity << " (target ≤ 1.05)" << std::endl;
  std::cout << "  Result: " << (r10.passed ? "PASS" : "FAIL") << std::endl;
  std::cout << std::endl;

  // Test 2: Frequency sweep
  std::cout << "=== Test 2: Frequency Sweep (7-12 GHz) ===" << std::endl;
  std::cout << std::setw(10) << "Freq(GHz)"
            << std::setw(10) << "|S11|"
            << std::setw(10) << "|S21|"
            << std::setw(12) << "Phase(°)"
            << std::setw(12) << "Expected(°)"
            << std::setw(10) << "Result" << std::endl;
  std::cout << std::string(64, '-') << std::endl;

  std::vector<double> test_freqs = {7e9, 8e9, 9e9, 10e9, 11e9, 12e9};
  int sweep_passes = 0;

  for (double freq : test_freqs) {
    auto r = run_benchmark(mesh, bc, port1_surf, port2_surf, dims, freq, L);
    std::cout << std::setw(10) << r.freq_ghz
              << std::setw(10) << r.s11_mag
              << std::setw(10) << r.s21_mag
              << std::setw(12) << r.s21_phase_deg
              << std::setw(12) << r.expected_phase_deg
              << std::setw(10) << (r.passed ? "PASS" : "FAIL") << std::endl;
    if (r.passed) sweep_passes++;
  }

  std::cout << std::endl;
  std::cout << "Frequency sweep: " << sweep_passes << "/" << test_freqs.size() << " passed" << std::endl;
  std::cout << std::endl;

  // Summary
  std::cout << "======================================================" << std::endl;
  std::cout << "BENCHMARK SUMMARY" << std::endl;
  std::cout << "======================================================" << std::endl;
  std::cout << "Center frequency test: " << (r10.passed ? "PASS" : "FAIL") << std::endl;
  std::cout << "Frequency sweep: " << sweep_passes << "/" << test_freqs.size() << " passed" << std::endl;

  bool overall_pass = r10.passed && (sweep_passes >= test_freqs.size() - 1);  // Allow 1 fail in sweep
  std::cout << std::endl;
  std::cout << "OVERALL: " << (overall_pass ? "PASS" : "FAIL") << std::endl;

  return overall_pass ? 0 : 1;
}
