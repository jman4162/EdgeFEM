// A PML section terminating a WR-90 guide must absorb the TE10 wave.
//
// examples/rect_waveguide.msh (z = 0 to 50 mm, Port1 at z = 0) has its last
// 20 mm retagged as a PML region; the far end (Port2 face, no BC: a PMC
// wall) is left open behind it. Port 1 is driven through the validated
// eigenmode path. Without a PML the PMC end reflects everything
// (|S11| = 1). An absorbing layer must give a passive, small |S11|; a layer
// with the wrong stretch sign is a gain medium and gives |S11| > 1.
//
// The isotropic stretch s = 1 - j sigma / omega makes eps = mu = s, so it is
// matched only for a plane wave: for TE10 the interface reflection is about
// alpha kc^2 / (2 beta^2) with alpha = sigma / omega (0.11 for alpha = 0.3
// at 10 GHz). The thresholds below allow that.

#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iostream>

#include "edgefem/maxwell.hpp"

using namespace edgefem;

namespace {
constexpr double c0 = 299792458.0;
constexpr int kPmlTag = 999;

double run(const MaxwellParams &p_in, const Mesh &mesh, const BC &bc,
           double freq) {
  MaxwellParams p = p_in;
  p.omega = 2.0 * M_PI * freq;
  RectWaveguidePort dims{0.02286, 0.01016};
  PortMode mode = solve_te10_mode(dims, freq);
  const double kc_sq = (M_PI / dims.a) * (M_PI / dims.a);
  WavePort wp1 = build_wave_port_2d(mesh, 2, mode, bc.dirichlet_edges, kc_sq);
  Eigen::MatrixXcd S = calculate_sparams_eigenmode(mesh, p, bc, {wp1});
  return std::abs(S(0, 0));
}
} // namespace

int main() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  const double z_pml = 0.030; // PML from z = 30 mm to the end (50 mm)
  int n_pml = 0;
  for (auto &tet : mesh.tets) {
    double zc = 0.0;
    for (int i = 0; i < 4; ++i)
      zc += mesh.nodes[mesh.nodeIndex.at(tet.conn[i])].xyz.z() / 4.0;
    if (zc > z_pml) {
      tet.phys = kPmlTag;
      ++n_pml;
    }
  }
  CHECK(n_pml > 0);

  const double freq = 10e9;
  const double omega = 2.0 * M_PI * freq;

  // Reference: no absorber, the PMC end reflects the whole wave.
  MaxwellParams p_none;
  const double s11_none = run(p_none, mesh, bc, freq);
  std::cout << "no PML:      |S11| = " << s11_none << "\n";
  CHECK(s11_none > 0.9);

  // Uniform isotropic stretch, alpha = sigma / omega = 0.3.
  MaxwellParams p_scalar;
  p_scalar.pml_regions.insert(kPmlTag);
  p_scalar.pml_sigma = 0.3 * omega;
  const double s11_scalar = run(p_scalar, mesh, bc, freq);
  std::cout << "scalar PML:  |S11| = " << s11_scalar << "\n";
  CHECK(s11_scalar <= 1.0);
  CHECK(s11_scalar < 0.25);

  // Graded tensor spec along z only (quadratic, 20 mm). The profile grows
  // toward the outer end of the domain; the inner face of the slab is the
  // interface with the guide and must start at sigma = 0.
  MaxwellParams p_tensor;
  PMLRegionSpec spec;
  spec.sigma_max = Eigen::Vector3d(0.0, 0.0, 1.0 * omega);
  spec.thickness = Eigen::Vector3d(0.0, 0.0, 0.020);
  spec.grading_order = 2.0;
  p_tensor.pml_tensor_regions.emplace(kPmlTag, spec);
  const double s11_tensor = run(p_tensor, mesh, bc, freq);
  std::cout << "graded PML:  |S11| = " << s11_tensor << "\n";
  CHECK(s11_tensor <= 1.0);
  CHECK(s11_tensor < 0.15);

  return edgefem_test::exit_code();
}
