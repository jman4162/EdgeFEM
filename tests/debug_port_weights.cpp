// Debug script to check port weight computation
#include <iostream>
#include <iomanip>
#include <cmath>
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

int main() {
  std::cout << std::setprecision(8);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;

  RectWaveguidePort dims{0.02286, 0.01016};
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortMode mode = solve_te10_mode(dims, freq);

  std::cout << "=== Mode Properties ===" << std::endl;
  std::cout << "kc = " << mode.kc << std::endl;
  std::cout << "omega*mu/kc² = " << omega * std::real(mode.mu) / (mode.kc * mode.kc) << std::endl;

  populate_te10_field(port1_surf, dims, mode);

  std::cout << "\n=== Mode Field (Hz) ===" << std::endl;
  std::cout << "First 5 values (should be real):" << std::endl;
  for (int i = 0; i < std::min(5, (int)mode.field.size()); ++i) {
    std::cout << "  Hz[" << i << "] = " << mode.field(i) << std::endl;
  }

  WavePort wp = build_wave_port(mesh, port1_surf, mode);

  std::cout << "\n=== Port Weights ===" << std::endl;
  std::cout << "Number of edges: " << wp.edges.size() << std::endl;

  // Check if weights are complex
  double real_sum = 0, imag_sum = 0;
  double real_sq_sum = 0, imag_sq_sum = 0;
  for (int i = 0; i < wp.weights.size(); ++i) {
    real_sum += std::real(wp.weights(i));
    imag_sum += std::imag(wp.weights(i));
    real_sq_sum += std::real(wp.weights(i)) * std::real(wp.weights(i));
    imag_sq_sum += std::imag(wp.weights(i)) * std::imag(wp.weights(i));
  }

  std::cout << "Sum of real parts: " << real_sum << std::endl;
  std::cout << "Sum of imag parts: " << imag_sum << std::endl;
  std::cout << "Sum of |real|²: " << real_sq_sum << std::endl;
  std::cout << "Sum of |imag|²: " << imag_sq_sum << std::endl;
  std::cout << "||w||² = " << real_sq_sum + imag_sq_sum << std::endl;

  std::cout << "\nFirst 10 weights:" << std::endl;
  for (int i = 0; i < std::min(10, (int)wp.weights.size()); ++i) {
    std::cout << "  w[" << i << "] = " << wp.weights(i) << std::endl;
  }

  // Expected: weights should be purely imaginary because:
  // Et = j * (ωμ/kc²) * ∇Hz
  // where ∇Hz is real, so Et is purely imaginary
  // and weight = ∫ Et · dl is imaginary

  std::cout << "\n=== Source Term Check ===" << std::endl;
  // Source = 2 * w / sqrt(Z0)
  std::complex<double> scale = 2.0 / std::sqrt(wp.mode.Z0);
  std::cout << "Source scale 2/sqrt(Z0) = " << scale << std::endl;

  std::complex<double> b_sum(0, 0);
  for (int i = 0; i < wp.weights.size(); ++i) {
    b_sum += scale * wp.weights(i);
  }
  std::cout << "Sum of source terms: " << b_sum << std::endl;
  std::cout << "Expected: purely imaginary if weights are imaginary" << std::endl;

  // Actually assemble and check RHS
  MaxwellParams p;
  p.omega = omega;
  std::vector<WavePort> ports{wp};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);

  double b_real_sum = 0, b_imag_sum = 0;
  double b_real_max = 0, b_imag_max = 0;
  for (int i = 0; i < asmbl.b.size(); ++i) {
    b_real_sum += std::real(asmbl.b(i));
    b_imag_sum += std::imag(asmbl.b(i));
    b_real_max = std::max(b_real_max, std::abs(std::real(asmbl.b(i))));
    b_imag_max = std::max(b_imag_max, std::abs(std::imag(asmbl.b(i))));
  }

  std::cout << "\n=== Assembled RHS ===" << std::endl;
  std::cout << "Sum of real parts: " << b_real_sum << std::endl;
  std::cout << "Sum of imag parts: " << b_imag_sum << std::endl;
  std::cout << "Max |real|: " << b_real_max << std::endl;
  std::cout << "Max |imag|: " << b_imag_max << std::endl;

  // Solve and check solution
  auto res = solve_linear(asmbl.A, asmbl.b, {});

  double x_real_max = 0, x_imag_max = 0;
  for (int i = 0; i < res.x.size(); ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      x_real_max = std::max(x_real_max, std::abs(std::real(res.x(i))));
      x_imag_max = std::max(x_imag_max, std::abs(std::imag(res.x(i))));
    }
  }

  std::cout << "\n=== Solution ===" << std::endl;
  std::cout << "Max |real|: " << x_real_max << std::endl;
  std::cout << "Max |imag|: " << x_imag_max << std::endl;
  std::cout << "Expected: solution should have both real and imag parts" << std::endl;

  // Port voltage
  std::complex<double> V1(0, 0);
  for (size_t k = 0; k < wp.edges.size(); ++k) {
    int e = wp.edges[k];
    if (!bc.dirichlet_edges.count(e)) {
      V1 += std::conj(wp.weights(k)) * res.x(e);
    }
  }

  std::cout << "\n=== Port Voltage ===" << std::endl;
  std::cout << "V1 = " << V1 << std::endl;
  std::cout << "V_inc = sqrt(Z0) = " << std::sqrt(wp.mode.Z0) << std::endl;
  std::cout << "V1 / V_inc = " << V1 / std::sqrt(wp.mode.Z0) << std::endl;

  return 0;
}
