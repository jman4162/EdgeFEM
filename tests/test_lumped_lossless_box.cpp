// Lumped port on a lossless closed box: |S11| must equal 1.
//
// Geometry: examples/parallel_plate_line.msh (air box 20 x 7.53 x 1 mm).
// Every boundary face is PEC except the x = 0 face, which carries the
// lumped port (E across the 1 mm gap, z direction). The port face is the
// only place energy can leave, the air is lossless (eps_r = 1) and the ABC
// is off.
//
// Reference: a lossless one-port has a purely reactive input impedance,
// Z = jX, so S11 = (jX - Z0)/(jX + Z0) has |S11| = 1 exactly, at any mesh
// density and any frequency that is not a box resonance. With PEC side
// walls, the lowest box resonance is near 20 GHz (c/2 * sqrt(1/w^2 +
// 1/(2L)^2)), well above the test frequencies 2, 5 and 8 GHz.
//
// Checked numerically with the gap-port readout from the 2026-10-12
// investigation (Z = j k0 eta0 w^T A0^-1 w / W^2, raw surface-integral
// weights w, strip width W = 7.53 mm): Z = j4.23, j10.87, j18.36 ohm and
// |S11| = 1 to 1e-12 at 2, 5, 8 GHz on this mesh.
//
// On commit 2c4d96d the port term w w^H / Z0 is real, so S11 is real with
// |S11| < 1; with normalize_port_weights the result is 0. Both are wrong.
//
// Three production paths are checked: calculate_sparams on the port as
// built, the documented normalize_port_weights + calculate_sparams recipe,
// and frequency_sweep.

#include "edgefem/bc.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/ports/lumped_port.hpp"
#include "edgefem/sweep.hpp"
#include "test_check.hpp"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

namespace {

constexpr int kTagPlates = 1;
constexpr int kTagPort1 = 2;
constexpr int kTagPort2Face = 3;
constexpr int kTagSides = 4;
constexpr double kTol = 1e-3;

BC pec_from_tags(const Mesh &mesh, const std::vector<int> &tags) {
  BC bc;
  for (int tag : tags) {
    BC t = build_edge_pec(mesh, tag);
    bc.dirichlet_edges.insert(t.dirichlet_edges.begin(),
                              t.dirichlet_edges.end());
  }
  return bc;
}

void report(const char *path, double f, std::complex<double> s11) {
  std::cout << std::setw(28) << std::left << path << " f = " << f / 1e9
            << " GHz  S11 = " << s11 << "  |S11| = " << std::abs(s11)
            << std::endl;
}

} // namespace

int main() {
  std::cout << std::setprecision(6);
  Mesh mesh = load_gmsh_v2("examples/parallel_plate_line.msh");
  BC bc = pec_from_tags(mesh, {kTagPlates, kTagPort2Face, kTagSides});

  LumpedPortConfig cfg;
  cfg.surface_tag = kTagPort1;
  cfg.z0 = 50.0;
  cfg.e_direction = Eigen::Vector3d(0, 0, 1);
  const WavePort port = build_lumped_port(mesh, cfg);

  SolveOptions opts;
  opts.use_direct = true;

  MaxwellParams base;
  base.use_abc = false; // closed box, no radiation boundary
  base.eps_r = 1.0;     // lossless air

  const std::vector<double> freqs = {2e9, 5e9, 8e9};

  for (double f : freqs) {
    MaxwellParams p = base;
    p.omega = 2.0 * M_PI * f;

    // Port as built by build_lumped_port.
    {
      std::vector<WavePort> ports{port};
      Eigen::MatrixXcd S = calculate_sparams(mesh, p, bc, ports, opts);
      report("calculate_sparams", f, S(0, 0));
      CHECK_NEAR(std::abs(S(0, 0)), 1.0, kTol);
    }

    // Documented recipe: normalize_port_weights, then calculate_sparams.
    {
      std::vector<WavePort> ports{port};
      normalize_port_weights(mesh, p, bc, ports, opts);
      Eigen::MatrixXcd S = calculate_sparams(mesh, p, bc, ports, opts);
      report("normalize+calculate_sparams", f, S(0, 0));
      CHECK_NEAR(std::abs(S(0, 0)), 1.0, kTol);
    }
  }

  // Frequency sweep over the same points.
  {
    std::vector<WavePort> ports{port};
    SweepResult sweep = frequency_sweep(mesh, base, bc, ports, freqs, opts);
    CHECK(sweep.S_matrices.size() == freqs.size());
    for (size_t i = 0; i < sweep.S_matrices.size(); ++i) {
      const std::complex<double> s11 = sweep.S_matrices[i](0, 0);
      report("frequency_sweep", freqs[i], s11);
      CHECK_NEAR(std::abs(s11), 1.0, kTol);
    }
  }

  return edgefem_test::exit_code();
}
