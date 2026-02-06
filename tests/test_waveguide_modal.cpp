#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

int main() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  const double freq = 10e9;
  MaxwellParams p;
  p.omega = 2 * M_PI * freq;

  PortSurfaceMesh port1 = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2 = extract_surface_mesh(mesh, 3);

  auto modes1 =
      solve_port_eigens(port1.mesh, 1, p.omega, 1.0, 1.0, ModePolarization::TE);
  auto modes2 =
      solve_port_eigens(port2.mesh, 1, p.omega, 1.0, 1.0, ModePolarization::TE);

  if (modes1.empty() || modes2.empty()) {
    std::cerr << "Failed to compute port modes" << std::endl;
    return 1;
  }

  WavePort wp1 = build_wave_port(mesh, port1, modes1.front());
  WavePort wp2 = build_wave_port(mesh, port2, modes2.front());
  std::vector<WavePort> ports{wp1, wp2};

  auto S = calculate_sparams(mesh, p, bc, ports);

  RectWaveguidePort dims{0.02286, 0.01016};
  double length = 0.05;
  auto analytic = straight_waveguide_sparams(dims, length, freq);

  assert(S.rows() == 2 && S.cols() == 2);
  assert(std::abs(S(0, 0)) < 0.2);
  assert(std::abs(S(1, 1)) < 0.2);
  assert(std::abs(std::abs(S(1, 0)) - std::abs(analytic.s21)) < 0.2);

  return 0;
}

