#include <cassert>
#include <cmath>
#include <iostream>

#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"

using namespace vectorem;

int main() {
  auto mesh = load_gmsh_v2("examples/wr90_port.msh");
  assert(!mesh.nodes.empty());
  assert(!mesh.tris.empty());
  assert(!mesh.boundary_lines.empty());

  // For TM modes with Ez=0 on boundary, lowest mode is TM11.
  const int num_modes = 1;
  double omega = 2 * M_PI * 10e9;
  auto modes =
      solve_port_eigens(mesh, num_modes, omega, 1.0, 1.0, ModePolarization::TE);

  if (modes.size() != 1) {
    std::cerr << "Expected 1 mode but found " << modes.size() << "\n";
    return 1;
  }
  const auto &mode = modes[0];
  std::cout << "Found mode with fc = " << mode.fc / 1e9 << " GHz\n";

  // Analytical TE10 cutoff frequency for rectangular waveguide
  const double a = 0.02286; // WR-90 width
  const double b = 0.01016; // WR-90 height
  const double c0 = 299792458.0;
  const double kc_te10 = M_PI / a;
  const double fc_te10 = kc_te10 * c0 / (2.0 * M_PI);

  std::cout << "Analytical TE10 fc = " << fc_te10 / 1e9 << " GHz\n";

  // Check if the computed fc is within 5% of the analytical value.
  assert(std::abs(mode.fc - fc_te10) / fc_te10 < 0.05);

  return 0;
}
