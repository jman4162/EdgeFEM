#include <iostream>
#include <vector>

#include "vectorem/io/touchstone.hpp"
#include "vectorem/maxwell.hpp"
#include "vectorem/mesh.hpp"

using namespace vectorem;

int main(int argc, char **argv) {
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  BC bc = build_edge_pec(mesh, 1);

  std::vector<LumpedPort> ports;
  ports.push_back({10, 50.0}); // Port 1 on edge 10
  ports.push_back({20, 50.0}); // Port 2 on edge 20

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
