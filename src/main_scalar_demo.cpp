#include <algorithm>
#include <iostream>
#include <string>

#include "vectorem/bc.hpp"
#include "vectorem/fem.hpp"
#include "vectorem/mesh.hpp"
#include "vectorem/solver.hpp"

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <mesh.msh> <pec_tag>\n";
    return 1;
  }
  std::string mesh_path = argv[1];
  int pec_tag = std::stoi(argv[2]);

  auto mesh = vectorem::load_gmsh_v2(mesh_path);
  auto bc = vectorem::build_scalar_pec(mesh, pec_tag);

  vectorem::ScalarHelmholtzParams p;
  p.k0 = 2.0 * M_PI * 1.0;
  p.eps_r = 1.0;

  auto asmbl = vectorem::assemble_scalar_helmholtz(mesh, p, bc);

  vectorem::VecC b = asmbl.b;
  if (!mesh.nodes.empty()) {
    int mid = static_cast<int>(mesh.nodes.size() / 2);
    b(mid) = std::complex<double>(1.0, 0.0);
  }

  vectorem::SolveOptions opt;
  opt.use_bicgstab = true;
  auto res = vectorem::solve_linear(asmbl.A, b, opt);

  std::cout << "Method: " << res.method << " iters=" << res.iters
            << " resid=" << res.residual << "\n";

  for (int i = 0; i < std::min<int>(5, res.x.size()); ++i) {
    std::cout << "u[" << i << "] = " << std::abs(res.x[i]) << "\n";
  }
  return 0;
}
