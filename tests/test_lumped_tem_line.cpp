// Two lumped ports at the ends of a matched TEM line: S11 = 0, S21 = 1.
//
// Geometry: examples/parallel_plate_line.msh, an air-filled section
// 0 <= x <= L, 0 <= y <= w, 0 <= z <= h with L = 20 mm, h = 1 mm,
// w = eta0 * h / 50 = 7.5346 mm.
//
// Boundary conditions: the plates z = 0 and z = h are PEC. The side walls
// y = 0 and y = w keep the natural boundary condition of the curl-curl
// weak form, n x H = 0 (PMC). Port 1 is the whole end face x = 0, port 2
// the whole end face x = L, both with E across the gap (z direction) and
// Z0 = 50 ohm.
//
// Reference (exact for the continuous problem):
//   Between PEC plates with PMC side walls the TEM field is uniform in y
//   and z: E = z_hat E0 exp(-j k0 x), H = -y_hat (E0 / eta0) exp(-j k0 x).
//   The gap voltage is V = E0 h and the plate current I = (E0 / eta0) w,
//   so Zc = V / I = eta0 h / w = 50 ohm by the choice of w.
//   Each port spans the full cross-section with a uniform sheet, so it
//   launches only this TEM field; no higher-order mode is excited.
//   A line of impedance Zc = Z0 and electrical length k0 L between two Z0
//   ports has (ABCD = [cos, j Zc sin; j sin / Zc, cos] of k0 L)
//     S11 = S22 = 0,   S21 = S12 = exp(-j k0 L),
//   with the exp(+j omega t) convention the code uses (its ABC adds
//   +j k0, and the eigenmode tests expect phase(S21) = -beta L).
//   Equivalent open-circuit Z-parameters: Z11 = -j Zc cot(k0 L),
//   Z21 = -j Zc / sin(k0 L).
//
// Checked numerically with the gap-port readout from the 2026-10-12
// investigation (Z_ij = j k0 eta0 w_i^T A0^-1 w_j / W^2 with raw
// surface-integral weights, W = w): at 2, 5, 8 GHz on this mesh,
// |S11| <= 4.9e-4, |S21| = 1.000000, |S12 - S21| ~ 1e-15, and
// phase(S21) within 0.11 deg of -k0 L. Z11 and Z21 agree with the
// formulas above to 0.1% (2 GHz) and 0.9% (8 GHz).
//
// The tolerances below leave margin over those numbers: |S11| < 0.02 and
// |S21| > 0.95, reciprocity 1e-6, losslessness |S11|^2 + |S21|^2 = 1
// within 1e-3, phase within 1 deg, and Z11, Z21 (converted from S with
// Z = Z0 (I + S)(I - S)^-1) within 2% of -j Zc cot(k0 L) and
// -j Zc / sin(k0 L). frequency_sweep must reproduce calculate_sparams to
// 1e-6 at the same frequencies.

#include "edgefem/bc.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/ports/lumped_port.hpp"
#include "edgefem/sweep.hpp"
#include "test_check.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

namespace {

constexpr double kC0 = 299792458.0;
constexpr double kLength = 20.0e-3;
constexpr double kZ0 = 50.0;
constexpr double kZc = 50.0; // eta0 h / w by the choice of w
constexpr int kTagPlates = 1;
constexpr int kTagPort1 = 2;
constexpr int kTagPort2 = 3;

WavePort make_port(const Mesh &mesh, int tag) {
  LumpedPortConfig cfg;
  cfg.surface_tag = tag;
  cfg.z0 = 50.0;
  cfg.e_direction = Eigen::Vector3d(0, 0, 1);
  return build_lumped_port(mesh, cfg);
}

void check_matched_line(const char *path, double f, const Eigen::MatrixXcd &S) {
  const double k0 = 2.0 * M_PI * f / kC0;
  const std::complex<double> s21_ref = std::polar(1.0, -k0 * kLength);
  const double phase_err_deg =
      std::abs(std::arg(S(1, 0) / s21_ref)) * 180.0 / M_PI;
  std::cout << path << " f = " << f / 1e9 << " GHz\n"
            << "  |S11| = " << std::abs(S(0, 0))
            << "  |S22| = " << std::abs(S(1, 1))
            << "  |S21| = " << std::abs(S(1, 0))
            << "  |S12| = " << std::abs(S(0, 1)) << "\n"
            << "  phase(S21) = " << std::arg(S(1, 0)) * 180.0 / M_PI
            << " deg, reference " << std::arg(s21_ref) * 180.0 / M_PI << " deg"
            << std::endl;

  CHECK(S.rows() == 2 && S.cols() == 2);
  CHECK(std::abs(S(0, 0)) < 0.02);
  CHECK(std::abs(S(1, 1)) < 0.02);
  CHECK(std::abs(S(1, 0)) > 0.95);
  CHECK(std::abs(S(0, 1)) > 0.95);
  CHECK_NEAR(S(0, 1), S(1, 0), 1e-6);
  CHECK_NEAR(std::norm(S(0, 0)) + std::norm(S(1, 0)), 1.0, 1e-3);
  CHECK_NEAR(std::norm(S(1, 1)) + std::norm(S(0, 1)), 1.0, 1e-3);
  CHECK(phase_err_deg < 1.0);

  // Open-circuit Z-parameters of the line section.
  const Eigen::MatrixXcd I2 = Eigen::MatrixXcd::Identity(2, 2);
  const Eigen::MatrixXcd Z = kZ0 * (I2 + S) * (I2 - S).inverse();
  const std::complex<double> j(0.0, 1.0);
  const std::complex<double> z11_ref = -j * kZc / std::tan(k0 * kLength);
  const std::complex<double> z21_ref = -j * kZc / std::sin(k0 * kLength);
  const double z11_rel = std::abs(Z(0, 0) - z11_ref) / std::abs(z11_ref);
  const double z21_rel = std::abs(Z(1, 0) - z21_ref) / std::abs(z21_ref);
  std::cout << "  Z11 = " << Z(0, 0) << " ohm, reference " << z11_ref
            << " (rel err " << z11_rel << ")\n"
            << "  Z21 = " << Z(1, 0) << " ohm, reference " << z21_ref
            << " (rel err " << z21_rel << ")" << std::endl;
  CHECK(z11_rel < 0.02);
  CHECK(z21_rel < 0.02);
}

} // namespace

int main() {
  std::cout << std::setprecision(6);
  Mesh mesh = load_gmsh_v2("examples/parallel_plate_line.msh");
  BC bc = build_edge_pec(mesh, kTagPlates); // side walls stay natural (PMC)

  const std::vector<WavePort> ports{make_port(mesh, kTagPort1),
                                    make_port(mesh, kTagPort2)};

  SolveOptions opts;
  opts.use_direct = true;

  MaxwellParams base;
  base.use_abc = false;
  base.eps_r = 1.0;

  const std::vector<double> freqs = {2e9, 5e9, 8e9};
  std::vector<Eigen::MatrixXcd> s_direct;

  for (double f : freqs) {
    MaxwellParams p = base;
    p.omega = 2.0 * M_PI * f;

    Eigen::MatrixXcd S = calculate_sparams(mesh, p, bc, ports, opts);
    check_matched_line("calculate_sparams", f, S);
    s_direct.push_back(S);

    std::vector<WavePort> normalized = ports;
    normalize_port_weights(mesh, p, bc, normalized, opts);
    Eigen::MatrixXcd Sn = calculate_sparams(mesh, p, bc, normalized, opts);
    check_matched_line("normalize+calculate_sparams", f, Sn);
  }

  SweepResult sweep = frequency_sweep(mesh, base, bc, ports, freqs, opts);
  CHECK(sweep.S_matrices.size() == freqs.size());
  for (size_t i = 0; i < sweep.S_matrices.size() && i < freqs.size(); ++i) {
    check_matched_line("frequency_sweep", freqs[i], sweep.S_matrices[i]);
    const double diff =
        (sweep.S_matrices[i] - s_direct[i]).cwiseAbs().maxCoeff();
    CHECK_NEAR(diff, 0.0, 1e-6);
  }

  return edgefem_test::exit_code();
}
