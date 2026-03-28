// ABC Scaling Factor Parametric Sweep
// Validates theoretical derivation that optimal coefficient is 0.5×jβ
//
// This test sweeps port_abc_scale from 0.1 to 2.0 and measures S-parameter
// accuracy to empirically verify the optimal value of 0.5.

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;

struct SweepResult {
  double alpha;       // ABC scaling factor
  double s11_mag;     // |S11|
  double s21_mag;     // |S21|
  double s21_phase;   // arg(S21) in degrees
  double phase_error; // Error vs analytical phase
  double passivity;   // |S11|^2 + |S21|^2
  double s21_error;   // |1 - |S21||
};

int main(int argc, char *argv[]) {
  std::cout << std::setprecision(6);
  std::cout << "=== ABC Scaling Factor Parametric Sweep ===" << std::endl;
  std::cout << "Purpose: Validate theoretical optimal value of alpha = 0.5"
            << std::endl
            << std::endl;

  // Load mesh and set up BCs
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  // Simulation parameters
  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double L = 0.05; // 50mm waveguide length

  // WR-90 waveguide dimensions
  RectWaveguidePort dims{0.02286, 0.01016};
  double kc = M_PI / dims.a;
  double kc_sq = kc * kc;
  double beta = std::sqrt(k0 * k0 - kc_sq);
  double Z0_te10 = omega * mu0 / beta;

  // Expected phase
  double expected_phase_rad = -beta * L;
  double expected_phase_deg = expected_phase_rad * 180.0 / M_PI;
  // Normalize to [-180, 180]
  while (expected_phase_deg < -180.0)
    expected_phase_deg += 360.0;
  while (expected_phase_deg > 180.0)
    expected_phase_deg -= 360.0;

  std::cout << "WR-90 Waveguide Parameters:" << std::endl;
  std::cout << "  Dimensions: " << dims.a * 1000 << " x " << dims.b * 1000
            << " mm" << std::endl;
  std::cout << "  Length: " << L * 1000 << " mm" << std::endl;
  std::cout << "  Frequency: " << freq / 1e9 << " GHz" << std::endl;
  std::cout << "  Beta: " << beta << " rad/m" << std::endl;
  std::cout << "  Z0: " << Z0_te10 << " ohms" << std::endl;
  std::cout << "  Expected phase: " << expected_phase_deg << " degrees"
            << std::endl
            << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Compute 3D FEM eigenvector for TE10 mode
  Eigen::VectorXd v_te10 =
      compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);

  // Create analytical mode parameters
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);

  // Build ports using eigenvector-based weights
  WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                  mode1, bc.dirichlet_edges);
  WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                  mode2, bc.dirichlet_edges);

  std::vector<WavePort> ports{wp1, wp2};

  // Sweep alpha values
  std::vector<double> alphas = {0.1,  0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
                                0.65, 0.7, 0.8, 0.9,  1.0, 1.2,  1.5, 2.0};

  std::vector<SweepResult> results;

  std::cout << "Sweeping ABC scaling factor alpha..." << std::endl;
  std::cout << std::setw(8) << "alpha" << std::setw(10) << "|S11|"
            << std::setw(10) << "|S21|" << std::setw(12) << "phase(S21)"
            << std::setw(12) << "phase_err" << std::setw(12) << "passivity"
            << std::setw(12) << "|S21|_err" << std::endl;
  std::cout << std::string(76, '-') << std::endl;

  for (double alpha : alphas) {
    MaxwellParams p;
    p.omega = omega;
    p.port_abc_scale = alpha;

    auto S = calculate_sparams_eigenmode(mesh, p, bc, ports);

    SweepResult r;
    r.alpha = alpha;
    r.s11_mag = std::abs(S(0, 0));
    r.s21_mag = std::abs(S(1, 0));
    r.s21_phase = std::arg(S(1, 0)) * 180.0 / M_PI;
    r.phase_error = std::abs(r.s21_phase - expected_phase_deg);
    if (r.phase_error > 180.0)
      r.phase_error = 360.0 - r.phase_error;
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    r.s21_error = std::abs(1.0 - r.s21_mag);

    results.push_back(r);

    std::cout << std::setw(8) << std::fixed << std::setprecision(2) << alpha
              << std::setw(10) << std::setprecision(4) << r.s11_mag
              << std::setw(10) << r.s21_mag << std::setw(12)
              << std::setprecision(1) << r.s21_phase << std::setw(12)
              << r.phase_error << std::setw(12) << std::setprecision(4)
              << r.passivity << std::setw(12) << std::setprecision(4)
              << r.s21_error << std::endl;
  }

  // Find optimal alpha (minimum |S21| error)
  double best_alpha = 0.0;
  double best_s21_error = 1e9;
  for (const auto &r : results) {
    if (r.s21_error < best_s21_error) {
      best_s21_error = r.s21_error;
      best_alpha = r.alpha;
    }
  }

  std::cout << std::endl;
  std::cout << "=== Results Summary ===" << std::endl;
  std::cout << "Optimal alpha (minimum |S21| error): " << best_alpha
            << std::endl;
  std::cout << "Best |S21| error: " << best_s21_error * 100 << "%" << std::endl;

  // Export to CSV for plotting
  std::string csv_path = "abc_scaling_sweep.csv";
  if (argc > 1) {
    csv_path = argv[1];
  }

  std::ofstream csv(csv_path);
  csv << "alpha,s11_mag,s21_mag,s21_phase,phase_error,passivity,s21_error"
      << std::endl;
  for (const auto &r : results) {
    csv << r.alpha << "," << r.s11_mag << "," << r.s21_mag << "," << r.s21_phase
        << "," << r.phase_error << "," << r.passivity << "," << r.s21_error
        << std::endl;
  }
  csv.close();
  std::cout << "Results exported to: " << csv_path << std::endl;

  // Verification: optimal should be near 0.5
  bool optimal_near_half = (best_alpha >= 0.4 && best_alpha <= 0.6);
  std::cout << std::endl;
  std::cout << "=== Verification ===" << std::endl;
  std::cout << "Theoretical prediction: alpha_opt = 0.5" << std::endl;
  std::cout << "Measured optimal: alpha = " << best_alpha << std::endl;
  std::cout << "Optimal in range [0.4, 0.6]: "
            << (optimal_near_half ? "PASS" : "FAIL") << std::endl;

  // Additional check: alpha=0.5 should give |S21| > 0.95
  double s21_at_half = 0.0;
  for (const auto &r : results) {
    if (std::abs(r.alpha - 0.5) < 0.01) {
      s21_at_half = r.s21_mag;
      break;
    }
  }
  bool good_transmission = s21_at_half > 0.95;
  std::cout << "|S21| at alpha=0.5: " << s21_at_half
            << (good_transmission ? " (PASS: > 0.95)" : " (FAIL: < 0.95)")
            << std::endl;

  return (optimal_near_half && good_transmission) ? 0 : 1;
}
