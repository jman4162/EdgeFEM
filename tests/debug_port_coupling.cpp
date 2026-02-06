// Debug why port weights don't couple to the propagating mode
// Key question: Does the port weight vector w align with the FEM eigenvector
// for TE10?

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;

  RectWaveguidePort dims{0.02286, 0.01016};

  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  std::cout << "=== Port Edge Analysis ===" << std::endl;
  std::cout << "Port 1: " << wp1.edges.size() << " edges" << std::endl;
  std::cout << "Port 2: " << wp2.edges.size() << " edges" << std::endl;

  // Count PEC vs free edges
  int pec_count1 = 0, free_count1 = 0;
  for (int e : wp1.edges) {
    if (bc.dirichlet_edges.count(e))
      pec_count1++;
    else
      free_count1++;
  }
  std::cout << "Port 1: " << pec_count1 << " PEC edges, " << free_count1
            << " free edges" << std::endl;

  int pec_count2 = 0, free_count2 = 0;
  for (int e : wp2.edges) {
    if (bc.dirichlet_edges.count(e))
      pec_count2++;
    else
      free_count2++;
  }
  std::cout << "Port 2: " << pec_count2 << " PEC edges, " << free_count2
            << " free edges" << std::endl;

  // Analyze port weight pattern
  std::cout << "\n=== Port Weight Pattern ===" << std::endl;
  double w1_norm_sq = wp1.weights.squaredNorm();
  double w2_norm_sq = wp2.weights.squaredNorm();
  std::cout << "||w1||² = " << w1_norm_sq << std::endl;
  std::cout << "||w2||² = " << w2_norm_sq << std::endl;

  // Find max weight magnitude
  double max_w1 = 0, max_w2 = 0;
  for (int i = 0; i < wp1.weights.size(); ++i) {
    max_w1 = std::max(max_w1, std::abs(wp1.weights(i)));
  }
  for (int i = 0; i < wp2.weights.size(); ++i) {
    max_w2 = std::max(max_w2, std::abs(wp2.weights(i)));
  }
  std::cout << "max|w1| = " << max_w1 << std::endl;
  std::cout << "max|w2| = " << max_w2 << std::endl;

  // Show first few weights (non-PEC only)
  std::cout << "\nPort 1 weights (first 10 free edges):" << std::endl;
  int shown = 0;
  for (int i = 0; i < (int)wp1.edges.size() && shown < 10; ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      std::cout << "  edge " << e << ": w = " << wp1.weights(i) << std::endl;
      shown++;
    }
  }

  // Solve the system and examine the solution
  std::cout << "\n=== FEM Solution Analysis ===" << std::endl;
  MaxwellParams p;
  p.omega = omega;

  std::vector<WavePort> ports{wp1, wp2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0); // Port 1 excited
  auto res = solve_linear(asmbl.A, asmbl.b, {});

  // Check solution norm
  double x_norm = res.x.norm();
  std::cout << "||x|| = " << x_norm << std::endl;

  // Check RHS norm
  double b_norm = asmbl.b.norm();
  std::cout << "||b|| = " << b_norm << std::endl;

  // Examine solution at port edges
  std::cout << "\nPort 1 solution (first 10 free edges):" << std::endl;
  shown = 0;
  double sum_wx1 = 0;
  for (int i = 0; i < (int)wp1.edges.size() && shown < 10; ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      std::complex<double> x_e = res.x(e);
      std::complex<double> w_e = wp1.weights(i);
      std::cout << "  edge " << e << ": x = " << x_e
                << ", w*x = " << std::conj(w_e) * x_e << std::endl;
      shown++;
    }
  }

  std::cout << "\nPort 2 solution (first 10 free edges):" << std::endl;
  shown = 0;
  for (int i = 0; i < (int)wp2.edges.size() && shown < 10; ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      std::complex<double> x_e = res.x(e);
      std::complex<double> w_e = wp2.weights(i);
      std::cout << "  edge " << e << ": x = " << x_e
                << ", w*x = " << std::conj(w_e) * x_e << std::endl;
      shown++;
    }
  }

  // Compute port voltages
  std::complex<double> V1 = 0, V2 = 0;
  for (int i = 0; i < (int)wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1 += std::conj(wp1.weights(i)) * res.x(e);
    }
  }
  for (int i = 0; i < (int)wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2 += std::conj(wp2.weights(i)) * res.x(e);
    }
  }

  std::complex<double> V_inc = std::sqrt(wp1.mode.Z0);
  std::complex<double> S11 = (V1 - V_inc) / V_inc;
  std::complex<double> S21 = V2 / V_inc;

  std::cout << "\n=== S-parameter Calculation ===" << std::endl;
  std::cout << "V_inc = " << V_inc << std::endl;
  std::cout << "V1 = " << V1 << std::endl;
  std::cout << "V2 = " << V2 << std::endl;
  std::cout << "S11 = " << S11 << " (|S11| = " << std::abs(S11) << ")"
            << std::endl;
  std::cout << "S21 = " << S21 << " (|S21| = " << std::abs(S21) << ")"
            << std::endl;

  // Expected values
  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::complex<double> S21_expected =
      std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\nExpected: S11 = 0, S21 = " << S21_expected << std::endl;

  // Ratio analysis
  std::cout << "\n=== Ratio Analysis ===" << std::endl;
  std::cout << "V1 / V_inc = " << V1 / V_inc << std::endl;
  std::cout << "V2 / V1 = " << V2 / V1 << std::endl;
  std::cout << "||w||² / Z0 = " << w1_norm_sq / std::real(wp1.mode.Z0)
            << std::endl;

  // Check the source term distribution
  std::cout << "\n=== Source Term Analysis ===" << std::endl;
  std::complex<double> scale = 2.0 / std::sqrt(wp1.mode.Z0);
  std::cout << "Source scale = 2/sqrt(Z0) = " << scale << std::endl;

  double b_port_norm = 0;
  for (int i = 0; i < (int)wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      std::complex<double> b_e = scale * wp1.weights(i);
      b_port_norm += std::norm(b_e);
    }
  }
  b_port_norm = std::sqrt(b_port_norm);
  std::cout << "||b_port|| = " << b_port_norm << std::endl;
  std::cout << "||b_total|| = " << b_norm << std::endl;
  std::cout << "Ratio ||b_port||/||b_total|| = " << b_port_norm / b_norm
            << std::endl;

  // Compare with matrix diagonal
  std::cout << "\n=== Matrix Diagonal at Port Edges ===" << std::endl;
  double avg_diag = 0;
  int diag_count = 0;
  for (int i = 0; i < (int)wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      double diag = std::abs(asmbl.A.coeff(e, e));
      avg_diag += diag;
      diag_count++;
      if (diag_count <= 5) {
        std::cout << "  A(" << e << "," << e << ") = " << asmbl.A.coeff(e, e)
                  << std::endl;
      }
    }
  }
  avg_diag /= diag_count;
  std::cout << "Average |A_ii| at port = " << avg_diag << std::endl;

  // The key ratio: source / diagonal determines solution magnitude
  std::cout << "\n=== Expected Solution Magnitude ===" << std::endl;
  std::cout << "If x ~ b/A_ii, then |x| ~ " << b_port_norm / avg_diag
            << std::endl;
  std::cout << "And V ~ ||w|| * |x| ~ "
            << std::sqrt(w1_norm_sq) * b_port_norm / avg_diag << std::endl;

  return 0;
}
