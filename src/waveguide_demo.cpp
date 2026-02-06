#include <iostream>
#include <vector>

#include "edgefem/io/touchstone.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"

using namespace edgefem;

int main(int argc, char **argv) {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  PortSurfaceMesh port1_surface = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surface = extract_surface_mesh(mesh, 3);

  double f0 = 8e9;  // Hz
  double f1 = 12e9; // Hz
  int N = 11;

  std::vector<double> freqs;
  std::vector<SParams2> s_params;

  std::cout << "Running frequency sweep...\n";
  for (int i = 0; i < N; ++i) {
    double f = f0 + (f1 - f0) * i / (N - 1);
    freqs.push_back(f);

    MaxwellParams p;
    p.omega = 2 * M_PI * f;

    auto modes1 =
        solve_port_eigens(port1_surface.mesh, 1, p.omega, 1.0, 1.0,
                          ModePolarization::TE);
    auto modes2 =
        solve_port_eigens(port2_surface.mesh, 1, p.omega, 1.0, 1.0,
                          ModePolarization::TE);
    std::vector<WavePort> ports;
    ports.push_back(build_wave_port(mesh, port1_surface, modes1.at(0)));
    ports.push_back(build_wave_port(mesh, port2_surface, modes2.at(0)));

    auto S = calculate_sparams(mesh, p, bc, ports);

    SParams2 s2;
    s2.s11 = S(0,0);
    s2.s21 = S(1,0);
    s2.s12 = S(0,1);
    s2.s22 = S(1,1);
    s_params.push_back(s2);

    std::cout << "f=" << f/1e9 << " GHz, S11=" << std::abs(s2.s11) << "\n";
  }

  write_touchstone("waveguide_sparams.s2p", freqs, s_params);
  std::cout << "Wrote waveguide_sparams.s2p\n";

  return 0;
}
