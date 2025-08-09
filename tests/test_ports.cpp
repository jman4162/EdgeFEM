#include <cassert>
#include <cmath>
#include <filesystem>

#include "vectorem/io/touchstone.hpp"
#include "vectorem/ports/port_eigensolve.hpp"

using namespace vectorem;

int main() {
  const RectWaveguidePort wg{0.02286, 0.01016}; // WR-90 dimensions in meters
  const double freq = 10e9;                     // 10 GHz
  const auto mode = solve_te10_mode(wg, freq);
  const double fc_ref = 299792458.0 / (2.0 * wg.a);
  assert(std::abs(mode.fc - fc_ref) < 1e6); // within 1 MHz
  const double eta0 = 376.730313668;
  const double z0_ref = eta0 / std::sqrt(1.0 - std::pow(fc_ref / freq, 2.0));
  assert(std::abs(mode.Z0 - z0_ref) / z0_ref < 0.01); // within 1%

  const double length = 0.1; // 10 cm straight section
  const auto s = straight_waveguide_sparams(wg, length, freq);
  assert(std::abs(s.s11) < 0.1); // |S11| < -20 dB
  assert(std::abs(std::abs(s.s21) - 1.0) < 1e-6);

  std::filesystem::create_directory("out");
  write_touchstone("out/test.s2p", {freq}, {s});
  assert(std::filesystem::exists("out/test.s2p"));

  return 0;
}
