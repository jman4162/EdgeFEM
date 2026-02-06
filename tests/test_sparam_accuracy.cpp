// Comprehensive S-parameter accuracy test
// Tests different normalization and ABC strategies to find the best approach
// Target: |S21| > 0.95, |S11| < 0.1 for matched 50mm WR-90 waveguide

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);
constexpr double eta0 = mu0 * c0;

struct TestResult {
  std::string name;
  std::complex<double> S11, S21;
  double passivity;
  double phase_error_deg;
};

void print_result(const TestResult &r, double expected_phase_deg) {
  double phase_deg = std::arg(r.S21) * 180.0 / M_PI;
  std::cout << std::setw(35) << std::left << r.name << " |S11|=" << std::setw(8)
            << std::abs(r.S11) << " |S21|=" << std::setw(8) << std::abs(r.S21)
            << " phase=" << std::setw(8) << phase_deg << "°"
            << " pass=" << std::setw(6) << r.passivity;
  if (std::abs(r.S21) > 0.5) {
    std::cout << " ✓";
  } else {
    std::cout << " ✗";
  }
  std::cout << std::endl;
}

int main() {
  std::cout << std::setprecision(4) << std::fixed;

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double L = 0.05; // 50mm waveguide

  RectWaveguidePort dims{0.02286, 0.01016}; // WR-90
  double kc_te10 = M_PI / dims.a;
  double beta = std::sqrt(k0 * k0 - kc_te10 * kc_te10);
  double expected_phase_deg = -beta * L * 180.0 / M_PI;

  std::cout << "=== WR-90 Waveguide S-Parameter Accuracy Test ===" << std::endl;
  std::cout << "Frequency: " << freq / 1e9 << " GHz" << std::endl;
  std::cout << "Length: " << L * 1000 << " mm" << std::endl;
  std::cout << "Expected: |S11| ≈ 0, |S21| ≈ 1, phase(S21) = "
            << expected_phase_deg << "°" << std::endl;
  std::cout << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Create mode parameters
  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  double Z0 = std::real(mode1.Z0);

  std::cout << "Mode Z0 = " << Z0 << " ohms" << std::endl;
  std::cout << "Mode beta = " << std::real(mode1.beta) << " rad/m" << std::endl;
  std::cout << std::endl;

  // Compute 3D FEM eigenvector
  double kc_te10_sq = kc_te10 * kc_te10;
  Eigen::VectorXd v_te10 =
      compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_te10_sq);

  // Build ports with different weight strategies
  std::vector<TestResult> results;

  // ========================================
  // Test 1: Analytical weights (baseline)
  // ========================================
  {
    populate_te10_field(port1_surf, dims, mode1);
    populate_te10_field(port2_surf, dims, mode2);
    WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
    WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

    // Normalize to ||w||² = sqrt(Z0)
    double target_norm_sq = std::sqrt(Z0);
    double scale1 = std::sqrt(target_norm_sq / wp1.weights.squaredNorm());
    double scale2 = std::sqrt(target_norm_sq / wp2.weights.squaredNorm());
    wp1.weights *= scale1;
    wp2.weights *= scale2;

    MaxwellParams p;
    p.omega = omega;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Analytical (||w||²=sqrt(Z0))";
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 2: Eigenvector weights (baseline)
  // ========================================
  {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    MaxwellParams p;
    p.omega = omega;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eigenvector (||w||²=sqrt(Z0))";
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 3: Eigenvector + normalize_port_weights
  // ========================================
  {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    MaxwellParams p;
    p.omega = omega;

    std::vector<WavePort> ports{wp1, wp2};

    std::cerr << "\n--- Normalizing port weights ---" << std::endl;
    normalize_port_weights(mesh, p, bc, ports);

    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eigenvector + norm (wAinvw=Z0/2)";
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 4-7: Different ABC types with eigenvector
  // ========================================
  std::vector<std::pair<PortABCType, std::string>> abc_types = {
      {PortABCType::Beta, "ABC: j*beta"},
      {PortABCType::BetaNorm, "ABC: j*beta/k0"},
      {PortABCType::ImpedanceMatch, "ABC: j*beta*sqrt(Z0/eta0)"},
      {PortABCType::ModalAdmittance, "ABC: j*omega*eps0/Z0"}};

  for (const auto &abc : abc_types) {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    MaxwellParams p;
    p.omega = omega;
    p.use_port_abc = true;
    p.port_abc_type = abc.first;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eigenvector + " + abc.second;
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 8-11: ABC types + normalization
  // ========================================
  for (const auto &abc : abc_types) {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    MaxwellParams p;
    p.omega = omega;
    p.use_port_abc = true;
    p.port_abc_type = abc.first;

    std::vector<WavePort> ports{wp1, wp2};
    normalize_port_weights(mesh, p, bc, ports);

    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Norm + " + abc.second;
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 12: Try different weight normalizations
  // ========================================
  // Test with ||w||² = Z0 instead of sqrt(Z0)
  {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    // Scale to ||w||² = Z0
    double target_norm_sq = Z0;
    double scale1 = std::sqrt(target_norm_sq / wp1.weights.squaredNorm());
    double scale2 = std::sqrt(target_norm_sq / wp2.weights.squaredNorm());
    wp1.weights *= scale1;
    wp2.weights *= scale2;

    MaxwellParams p;
    p.omega = omega;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eigenvector (||w||²=Z0)";
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // Test with ||w||² = 1
  {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    // Scale to ||w||² = 1
    double scale1 = 1.0 / wp1.weights.norm();
    double scale2 = 1.0 / wp2.weights.norm();
    wp1.weights *= scale1;
    wp2.weights *= scale2;

    MaxwellParams p;
    p.omega = omega;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eigenvector (||w||²=1)";
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // ========================================
  // Test 14-17: Different power scales with ABC
  // ========================================
  std::vector<std::pair<double, std::string>> power_scales = {
      {1.0, "||w||²=1"},
      {Z0, "||w||²=Z0"},
      {std::sqrt(Z0), "||w||²=sqrt(Z0)"},
      {Z0 / 2.0, "||w||²=Z0/2"}};

  for (const auto &ps : power_scales) {
    WavePort wp1 = build_wave_port_from_eigenvector(mesh, port1_surf, v_te10,
                                                    mode1, bc.dirichlet_edges);
    WavePort wp2 = build_wave_port_from_eigenvector(mesh, port2_surf, v_te10,
                                                    mode2, bc.dirichlet_edges);

    // Scale to target ||w||²
    double target = ps.first;
    double scale1 = std::sqrt(target / wp1.weights.squaredNorm());
    double scale2 = std::sqrt(target / wp2.weights.squaredNorm());
    wp1.weights *= scale1;
    wp2.weights *= scale2;

    MaxwellParams p;
    p.omega = omega;
    p.use_port_abc = true;
    p.port_abc_type = PortABCType::Beta;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    TestResult r;
    r.name = "Eig+ABC " + ps.second;
    r.S11 = S(0, 0);
    r.S21 = S(1, 0);
    r.passivity = std::norm(S(0, 0)) + std::norm(S(1, 0));
    results.push_back(r);
  }

  // Print results
  std::cout << "=== Results ===" << std::endl;
  std::cout << std::string(90, '-') << std::endl;
  for (const auto &r : results) {
    print_result(r, expected_phase_deg);
  }
  std::cout << std::string(90, '-') << std::endl;

  // Summary
  std::cout << "\nTarget: |S11| < 0.1, |S21| > 0.95, passivity ≤ 1.01"
            << std::endl;

  bool any_pass = false;
  for (const auto &r : results) {
    if (std::abs(r.S21) > 0.95 && std::abs(r.S11) < 0.1 &&
        r.passivity <= 1.01) {
      std::cout << "PASS: " << r.name << std::endl;
      any_pass = true;
    }
  }

  if (!any_pass) {
    std::cout << "No configuration meets all criteria yet." << std::endl;

    // Find best transmission
    double best_s21 = 0;
    std::string best_name;
    for (const auto &r : results) {
      if (std::abs(r.S21) > best_s21) {
        best_s21 = std::abs(r.S21);
        best_name = r.name;
      }
    }
    std::cout << "Best |S21| = " << best_s21 << " from: " << best_name
              << std::endl;
  }

  return any_pass ? 0 : 1;
}
