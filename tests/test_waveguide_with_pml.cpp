// Test waveguide with PML regions at the ends to support traveling waves
// Instead of wave ports with matched loads, use PML to absorb outgoing waves
// Then use a source in the interior to excite the TE10 mode

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  // For this test, we need a mesh with PML regions
  // The current rect_waveguide.msh doesn't have PML regions
  // Let's instead try with the existing mesh and add artificial damping at
  // ports

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

  std::cout << "=== Testing Alternative Formulations ===" << std::endl;
  std::cout << "Operating frequency: 10 GHz" << std::endl;
  std::cout << "Expected TE10 beta = " << std::real(mode1.beta) << " rad/m"
            << std::endl;

  MaxwellParams p;
  p.omega = omega;

  // Test 1: Standard wave port formulation (baseline)
  std::cout << "\n=== Test 1: Standard Wave Port ===" << std::endl;
  std::vector<WavePort> ports{wp1, wp2};
  auto S1 = calculate_sparams(mesh, p, bc, ports);
  std::cout << "S11 = " << S1(0, 0) << " (|S11| = " << std::abs(S1(0, 0)) << ")"
            << std::endl;
  std::cout << "S21 = " << S1(1, 0) << " (|S21| = " << std::abs(S1(1, 0)) << ")"
            << std::endl;

  // Test 2: Add imaginary damping to port edges (like ABC but stronger)
  std::cout << "\n=== Test 2: Port with ABC Damping ===" << std::endl;

  // Manually assemble with modified port treatment
  auto asmbl2 = assemble_maxwell(mesh, p, bc, ports, 0);

  // Add additional imaginary damping to port edges
  // For proper absorption, damping should be ~ jωμ/Z0 × mode_amplitude
  double Z0 = std::real(mode1.Z0);
  double beta = std::real(mode1.beta);
  std::complex<double> j(0, 1);

  // Damping coefficient for ABC: jk where k = omega/c
  // For waveguide mode, use jβ instead
  std::complex<double> damping = j * beta;

  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      asmbl2.A.coeffRef(e, e) += damping;
    }
  }
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      asmbl2.A.coeffRef(e, e) += damping;
    }
  }

  auto res2 = solve_linear(asmbl2.A, asmbl2.b, {});

  // Extract S-parameters manually
  std::complex<double> V1_2 = 0, V2_2 = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1_2 += std::conj(wp1.weights(i)) * res2.x(e);
    }
  }
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2_2 += std::conj(wp2.weights(i)) * res2.x(e);
    }
  }

  std::complex<double> V_inc = std::sqrt(Z0);
  std::complex<double> S11_2 = (V1_2 - V_inc) / V_inc;
  std::complex<double> S21_2 = V2_2 / V_inc;

  std::cout << "S11 = " << S11_2 << " (|S11| = " << std::abs(S11_2) << ")"
            << std::endl;
  std::cout << "S21 = " << S21_2 << " (|S21| = " << std::abs(S21_2) << ")"
            << std::endl;

  // Test 3: Use mass-matrix-based damping (more physical ABC)
  std::cout << "\n=== Test 3: Mass-weighted ABC ===" << std::endl;

  // For a proper first-order ABC, the damping should be proportional to the
  // mass matrix Add jβ × M_port to the system matrix

  auto asmbl3 = assemble_maxwell(mesh, p, bc, ports, 0);

  // For simplicity, use diagonal mass approximation
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      // Estimate local mass (proportional to edge length squared / volume)
      const auto &edge = mesh.edges[e];
      const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
      double L_edge = (p1 - p0).norm();
      double M_local = L_edge * L_edge * 0.001; // Rough estimate

      asmbl3.A.coeffRef(e, e) += damping * M_local;
    }
  }
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      const auto &edge = mesh.edges[e];
      const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
      double L_edge = (p1 - p0).norm();
      double M_local = L_edge * L_edge * 0.001;

      asmbl3.A.coeffRef(e, e) += damping * M_local;
    }
  }

  auto res3 = solve_linear(asmbl3.A, asmbl3.b, {});

  std::complex<double> V1_3 = 0, V2_3 = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1_3 += std::conj(wp1.weights(i)) * res3.x(e);
    }
  }
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2_3 += std::conj(wp2.weights(i)) * res3.x(e);
    }
  }

  std::complex<double> S11_3 = (V1_3 - V_inc) / V_inc;
  std::complex<double> S21_3 = V2_3 / V_inc;

  std::cout << "S11 = " << S11_3 << " (|S11| = " << std::abs(S11_3) << ")"
            << std::endl;
  std::cout << "S21 = " << S21_3 << " (|S21| = " << std::abs(S21_3) << ")"
            << std::endl;

  // Test 4: Try scaling the ABC damping
  std::cout << "\n=== Test 4: Scaled ABC Damping (sweep) ===" << std::endl;

  for (double scale : {0.1, 1.0, 10.0, 100.0, 1000.0}) {
    auto asmbl4 = assemble_maxwell(mesh, p, bc, ports, 0);

    std::complex<double> scaled_damping = damping * scale;

    for (size_t i = 0; i < wp1.edges.size(); ++i) {
      int e = wp1.edges[i];
      if (!bc.dirichlet_edges.count(e)) {
        asmbl4.A.coeffRef(e, e) += scaled_damping;
      }
    }
    for (size_t i = 0; i < wp2.edges.size(); ++i) {
      int e = wp2.edges[i];
      if (!bc.dirichlet_edges.count(e)) {
        asmbl4.A.coeffRef(e, e) += scaled_damping;
      }
    }

    auto res4 = solve_linear(asmbl4.A, asmbl4.b, {});

    std::complex<double> V1_4 = 0, V2_4 = 0;
    for (size_t i = 0; i < wp1.edges.size(); ++i) {
      int e = wp1.edges[i];
      if (!bc.dirichlet_edges.count(e)) {
        V1_4 += std::conj(wp1.weights(i)) * res4.x(e);
      }
    }
    for (size_t i = 0; i < wp2.edges.size(); ++i) {
      int e = wp2.edges[i];
      if (!bc.dirichlet_edges.count(e)) {
        V2_4 += std::conj(wp2.weights(i)) * res4.x(e);
      }
    }

    std::complex<double> S11_4 = (V1_4 - V_inc) / V_inc;
    std::complex<double> S21_4 = V2_4 / V_inc;

    double passivity = std::norm(S11_4) + std::norm(S21_4);

    std::cout << "Scale " << scale << ": |S11|=" << std::abs(S11_4)
              << ", |S21|=" << std::abs(S21_4) << ", passivity=" << passivity
              << std::endl;
  }

  // Expected values
  double L = 0.05;
  std::complex<double> S21_expected =
      std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\nExpected: S11 = 0, S21 = " << S21_expected << std::endl;

  return 0;
}
