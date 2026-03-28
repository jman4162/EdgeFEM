#include "edgefem/bc.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/ports/lumped_port.hpp"
#include "edgefem/solver.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace edgefem;

int main() {
  // Load a simple mesh that has surface tags we can use
  // Use the cube_cavity mesh which has known surface tags
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");

  std::cout << "Mesh loaded: " << mesh.tets.size() << " tets, "
            << mesh.tris.size() << " tris, " << mesh.edges.size() << " edges\n";

  // Find available surface tags
  auto tags = list_physical_tags(mesh);
  std::cout << "Surface tags:";
  for (int t : tags.surface_tags)
    std::cout << " " << t;
  std::cout << "\n";

  if (tags.surface_tags.empty()) {
    std::cout << "SKIP: no surface tags in mesh\n";
    return 0;
  }

  // Use the first available surface tag
  int port_tag = *tags.surface_tags.begin();
  std::cout << "Using surface tag " << port_tag << " for lumped port\n";

  // Test 1: Build lumped port
  LumpedPortConfig config;
  config.surface_tag = port_tag;
  config.z0 = 50.0;
  config.e_direction = Eigen::Vector3d(0, 0, 1);

  WavePort port = build_lumped_port(mesh, config);

  std::cout << "Port edges: " << port.edges.size() << "\n";
  std::cout << "Port Z0: " << port.mode.Z0 << "\n";

  assert(port.edges.size() > 0 && "Lumped port must have edges");
  assert(port.weights.size() == static_cast<int>(port.edges.size()) &&
         "Weight count must match edge count");
  assert(std::abs(port.mode.Z0 - 50.0) < 1e-10 && "Z0 must be 50 ohms");
  assert(port.mode.beta == 0.0 && "Lumped port beta must be 0");

  // Test 2: Weights are non-zero
  double w_norm = port.weights.norm();
  std::cout << "Weight norm: " << w_norm << "\n";
  assert(w_norm > 1e-10 && "Weights must be non-zero");

  // Test 3: Assembly with lumped port doesn't produce singular system
  BC bc;
  // Add PEC on a different surface if available
  for (int t : tags.surface_tags) {
    if (t != port_tag) {
      BC bc_t = build_edge_pec(mesh, t);
      for (int e : bc_t.dirichlet_edges)
        bc.dirichlet_edges.insert(e);
      break;
    }
  }

  MaxwellParams params;
  params.omega = 2 * M_PI * 10e9; // 10 GHz

  std::vector<WavePort> ports = {port};

  auto assembly = assemble_maxwell(mesh, params, bc, ports, 0);

  std::cout << "System size: " << assembly.A.rows() << " x "
            << assembly.A.cols() << "\n";
  std::cout << "RHS norm: " << assembly.b.norm() << "\n";

  assert(assembly.A.rows() > 0 && "System matrix must be non-empty");
  assert(assembly.b.norm() > 1e-15 &&
         "RHS must be non-zero when port is active");

  // Test 4: System is solvable
  auto result = solve_linear(assembly.A, assembly.b);
  VecC solution = result.x;
  double sol_norm = solution.norm();
  std::cout << "Solution norm: " << sol_norm << "\n";
  assert(sol_norm > 1e-20 && "Solution must be non-zero");
  assert(std::isfinite(sol_norm) && "Solution must be finite");

  std::cout << "\nAll lumped port tests passed!\n";
  return 0;
}
