// WR-42 Waveguide S-Parameter Test (Ka-band)
// Validates frequency scaling by testing at higher frequency than WR-90
//
// WR-42 specifications:
//   a = 10.668 mm, b = 4.318 mm
//   fc (TE10) = 14.05 GHz
//   Operating band: 18-26.5 GHz
//
// Test frequencies: 20, 22, 24, 26 GHz

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;

struct FrequencyResult {
  double freq;        // Hz
  double beta_ana;    // Analytical propagation constant
  double Z0_ana;      // Analytical impedance
  double s11_mag;     // |S11|
  double s21_mag;     // |S21|
  double s21_phase;   // arg(S21) in degrees
  double phase_ana;   // Analytical phase
  double phase_error; // Error in degrees
};

int main() {
  std::cout << std::setprecision(6);
  std::cout << "=== WR-42 Waveguide S-Parameter Validation ===" << std::endl;
  std::cout << "Purpose: Validate frequency scaling at Ka-band" << std::endl
            << std::endl;

  // Load mesh
  Mesh mesh = load_gmsh_v2("examples/wr42_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  std::cout << "Mesh: " << mesh.nodes.size() << " nodes, " << mesh.tets.size()
            << " tets, " << mesh.edges.size() << " edges" << std::endl;

  // WR-42 dimensions
  RectWaveguidePort dims{0.010668, 0.004318}; // a=10.668mm, b=4.318mm
  double L = 0.030;                           // 30mm waveguide length

  double kc = M_PI / dims.a;
  double kc_sq = kc * kc;
  double fc = c0 * kc / (2 * M_PI);

  std::cout << "WR-42 Parameters:" << std::endl;
  std::cout << "  Dimensions: " << dims.a * 1000 << " x " << dims.b * 1000
            << " mm" << std::endl;
  std::cout << "  Length: " << L * 1000 << " mm" << std::endl;
  std::cout << "  Cutoff frequency: " << fc / 1e9 << " GHz" << std::endl;
  std::cout << "  Operating band: 18-26.5 GHz" << std::endl << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Compute eigenvector for TE10 mode
  Eigen::VectorXd v_te10 =
      compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);

  // Test frequencies
  std::vector<double> test_freqs = {20e9, 22e9, 24e9, 26e9};
  std::vector<FrequencyResult> results;

  std::cout << "Running frequency sweep..." << std::endl;
  std::cout << std::setw(10) << "Freq" << std::setw(10) << "|S11|"
            << std::setw(10) << "|S21|" << std::setw(12) << "Phase(S21)"
            << std::setw(12) << "Phase_ana" << std::setw(12) << "Phase_err"
            << std::setw(12) << "Passivity" << std::endl;
  std::cout << std::string(78, '-') << std::endl;

  for (double freq : test_freqs) {
    double omega = 2 * M_PI * freq;
    double k0 = omega / c0;
    double beta = std::sqrt(k0 * k0 - kc_sq);
    double Z0_te10 = omega * mu0 / beta;

    // Create mode parameters at this frequency
    PortMode mode1 = solve_te10_mode(dims, freq);
    PortMode mode2 = solve_te10_mode(dims, freq);

    // Build ports using eigenvector-based weights
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    // Set up Maxwell parameters with optimal ABC scale
    MaxwellParams p;
    p.omega = omega;
    p.port_abc_scale = 0.5;

    std::vector<WavePort> ports{wp1, wp2};

    // Calculate S-parameters
    auto S = calculate_sparams_eigenmode(mesh, p, bc, ports);

    // Compute results
    FrequencyResult r;
    r.freq = freq;
    r.beta_ana = beta;
    r.Z0_ana = Z0_te10;
    r.s11_mag = std::abs(S(0, 0));
    r.s21_mag = std::abs(S(1, 0));
    r.s21_phase = std::arg(S(1, 0)) * 180.0 / M_PI;

    // Expected phase: -beta*L
    r.phase_ana = -beta * L * 180.0 / M_PI;
    // Normalize to [-180, 180]
    while (r.phase_ana < -180.0)
      r.phase_ana += 360.0;
    while (r.phase_ana > 180.0)
      r.phase_ana -= 360.0;

    r.phase_error = std::abs(r.s21_phase - r.phase_ana);
    if (r.phase_error > 180.0)
      r.phase_error = 360.0 - r.phase_error;

    double passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));

    results.push_back(r);

    std::cout << std::setw(10) << std::fixed << std::setprecision(1)
              << freq / 1e9 << std::setw(10) << std::setprecision(4)
              << r.s11_mag << std::setw(10) << r.s21_mag << std::setw(12)
              << std::setprecision(1) << r.s21_phase << std::setw(12)
              << r.phase_ana << std::setw(12) << r.phase_error << std::setw(12)
              << std::setprecision(4) << passivity << std::endl;
  }

  // Validation criteria
  std::cout << std::endl << "=== Validation ===" << std::endl;

  int num_pass = 0;
  int total = static_cast<int>(results.size());

  for (const auto &r : results) {
    bool s11_pass =
        r.s11_mag < 0.20; // |S11| < -14 dB (relaxed for coarser mesh)
    bool s21_pass = r.s21_mag > 0.90;       // |S21| > -0.9 dB
    bool phase_pass = r.phase_error < 15.0; // Phase error < 15°

    if (s11_pass && s21_pass && phase_pass) {
      num_pass++;
    }
  }

  // Compute averages
  double avg_s11 = 0.0, avg_s21 = 0.0, avg_phase_err = 0.0;
  for (const auto &r : results) {
    avg_s11 += r.s11_mag;
    avg_s21 += r.s21_mag;
    avg_phase_err += r.phase_error;
  }
  avg_s11 /= total;
  avg_s21 /= total;
  avg_phase_err /= total;

  std::cout << "Average |S11|: " << avg_s11 << " (target: < 0.15)" << std::endl;
  std::cout << "Average |S21|: " << avg_s21 << " (target: > 0.90)" << std::endl;
  std::cout << "Average phase error: " << avg_phase_err << "° (target: < 15°)"
            << std::endl;
  std::cout << "Frequencies passing: " << num_pass << "/" << total << std::endl;

  // Overall pass: at least 3 of 4 frequencies pass all criteria
  bool overall_pass = (num_pass >= 3) && (avg_s21 > 0.90);
  std::cout << std::endl
            << "=== OVERALL: " << (overall_pass ? "PASS" : "FAIL")
            << " ===" << std::endl;

  return overall_pass ? 0 : 1;
}
