// Test S-parameters at different frequencies to understand FEM behavior
// Lower frequencies should work better if the issue is K >> k0²M
//
// The ports are modal (build_wave_port with analytic TE10 weights), so
// calculate_sparams and frequency_sweep use the 1.0.0 modal port terms,
// whose S-parameters are not validated (see test_modal_port_1p0).
//
// Pass/fail, default run: frequency_sweep (K and M assembled once, port
// terms added per frequency) returns the same S-matrix as calculate_sparams
// with the same ports at every frequency, to 1e-9.
//
// With the argument "passivity" (ctest entry test_frequency_sweep_passivity,
// label known_issue): every column of S must satisfy
// |S11|^2 + |S21|^2 <= 1.01, as a lossless guide requires. The 1.0.0 modal
// terms fail this: on 2026-10-15 the sums were 12.3 at 7 GHz, 5.80 at
// 8 GHz, 2.65 at 10 GHz and 10.3 at 20 GHz (|S11| = 3.51, 2.41, 1.63,
// 2.71). Above 13.1 GHz (TE20 cutoff) the guide is multimode, so the
// single-mode match |S11| ~ 0, |S21| ~ 1 is not expected there; passivity
// still is.

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include "edgefem/sweep.hpp"
#include "test_check.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main(int argc, char **argv) {
  const bool check_passivity = argc > 1 && std::string(argv[1]) == "passivity";
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  RectWaveguidePort dims{0.02286, 0.01016};
  double fc = c0 / (2 * dims.a); // TE10 cutoff ~ 6.56 GHz

  std::cout << "TE10 cutoff frequency: " << fc / 1e9 << " GHz" << std::endl;
  std::cout << "\nSweeping frequency from 7 GHz to 30 GHz:" << std::endl;
  std::cout
      << "freq(GHz)  k0²      K/(k0²M)  |S11|    |S21|    Phase(S21)  Status"
      << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  std::vector<double> freqs = {7e9, 8e9, 10e9, 15e9, 20e9, 25e9, 30e9};

  // Direct solves: the iterative default took about 3 minutes here.
  SolveOptions opts;
  opts.use_direct = true;

  for (double freq : freqs) {
    if (freq <= fc) {
      std::cout << freq / 1e9 << " GHz: Below cutoff, skipping" << std::endl;
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
      if (bc.dirichlet_edges.count(i))
        continue;
      K_avg += std::real(asmbl.A.coeff(i, i));
      count++;
    }
    K_avg /= count;

    // Estimate M from A = K - k0²M => M = (K - A)/k0² ≈ K/k0² (since A ≈ K for
    // low freq) Actually just use the formula: ratio = K_avg / (some estimated
    // k0²M) For simplicity, use K_avg directly as a proxy for the K
    // contribution
    double km_ratio = K_avg / k0_sq; // Rough estimate

    // Calculate S-parameters
    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports, opts);

    double s11_mag = std::abs(S(0, 0));
    double s21_mag = std::abs(S(1, 0));
    double s21_phase = std::arg(S(1, 0)) * 180 / M_PI;

    // Expected S21 phase = -beta * L * 180/pi
    double L = 0.05;
    double s21_phase_expected = -beta * L * 180 / M_PI;
    while (s21_phase_expected < -180)
      s21_phase_expected += 360;

    // Check passivity
    double power_sum = s11_mag * s11_mag + s21_mag * s21_mag;
    std::string status = (power_sum <= 1.01) ? "PASS" : "FAIL";
    if (check_passivity)
      CHECK(power_sum <= 1.01);

    std::cout << std::setw(5) << freq / 1e9 << "      " << std::setw(8) << k0_sq
              << "  " << std::setw(8) << km_ratio << "  " << std::setw(7)
              << s11_mag << "  " << std::setw(7) << s21_mag << "  "
              << std::setw(10) << s21_phase << "°  " << status << std::endl;
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
    auto S = calculate_sparams(mesh, p, bc, ports, opts);

    double s11_mag = std::abs(S(0, 0));
    double s21_mag = std::abs(S(1, 0));

    std::cout << freq / 1e9 << " GHz: |S11|=" << s11_mag
              << ", |S21|=" << s21_mag
              << ", |S11|²+|S21|²=" << s11_mag * s11_mag + s21_mag * s21_mag
              << std::endl;
    if (check_passivity)
      CHECK(s11_mag * s11_mag + s21_mag * s21_mag <= 1.01);
  }

  // Same ports (built at 10 GHz) through frequency_sweep, which assembles
  // K and M once and adds the port terms per frequency.
  std::cout << "\n=== frequency_sweep (ports built at 10 GHz) ===" << std::endl;
  {
    const double f_port = 10e9;
    PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
    PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);
    PortMode mode1 = solve_te10_mode(dims, f_port);
    PortMode mode2 = solve_te10_mode(dims, f_port);
    populate_te10_field(port1_surf, dims, mode1);
    populate_te10_field(port2_surf, dims, mode2);
    std::vector<WavePort> ports{build_wave_port(mesh, port1_surf, mode1),
                                build_wave_port(mesh, port2_surf, mode2)};

    MaxwellParams p;
    p.omega = 2 * M_PI * f_port;
    SweepResult sweep = frequency_sweep(mesh, p, bc, ports, freqs, opts);
    CHECK(sweep.S_matrices.size() == freqs.size());
    for (size_t i = 0; i < sweep.S_matrices.size(); ++i) {
      const Eigen::MatrixXcd &S = sweep.S_matrices[i];
      const double col1 = std::norm(S(0, 0)) + std::norm(S(1, 0));
      const double col2 = std::norm(S(1, 1)) + std::norm(S(0, 1));
      std::cout << freqs[i] / 1e9 << " GHz: |S11|=" << std::abs(S(0, 0))
                << ", |S21|=" << std::abs(S(1, 0)) << ", |S11|²+|S21|²=" << col1
                << ", |S22|²+|S12|²=" << col2 << std::endl;
      if (check_passivity) {
        CHECK(col1 <= 1.01);
        CHECK(col2 <= 1.01);
      }

      // Same ports through calculate_sparams at this frequency.
      MaxwellParams pf;
      pf.omega = 2 * M_PI * freqs[i];
      const Eigen::MatrixXcd S_ref =
          calculate_sparams(mesh, pf, bc, ports, opts);
      const double scale = std::max(1.0, S_ref.cwiseAbs().maxCoeff());
      const double diff = (S - S_ref).cwiseAbs().maxCoeff();
      std::cout << "  max |S_sweep - S_calculate_sparams| = " << diff
                << std::endl;
      CHECK(diff <= 1e-9 * scale);
    }
  }

  return edgefem_test::exit_code();
}
