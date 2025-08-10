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
  auto modes = solve_port_eigens(mesh, num_modes);

  assert(modes.size() == 1);
  const auto& mode = modes[0];
  std::cout << "Found mode with fc = " << mode.fc / 1e9 << " GHz\n";

  // Analytical TM11 cutoff frequency for rectangular waveguide
  const double a = 0.02286; // WR-90 width
  const double b = 0.01016; // WR-90 height
  const double c0 = 299792458.0;
  const double kc_tm11 = std::sqrt(std::pow(M_PI / a, 2) + std::pow(M_PI / b, 2));
  const double fc_tm11 = kc_tm11 * c0 / (2.0 * M_PI);

  std::cout << "Analytical TM11 fc = " << fc_tm11 / 1e9 << " GHz\n";

  // Check if the computed fc is within 5% of the analytical value.
  // The mesh is coarse, so we don't expect super high accuracy.
  assert(std::abs(mode.fc - fc_tm11) / fc_tm11 < 0.05);

  return 0;
}
