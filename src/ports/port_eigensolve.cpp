#include "vectorem/ports/port_eigensolve.hpp"

#include <cmath>

namespace vectorem {

namespace {
constexpr double c0 = 299792458.0;       // speed of light in vacuum (m/s)
constexpr double eta0 = 376.730313668;   // impedance of free space (ohms)
}

PortMode solve_te10_mode(const RectWaveguidePort &port, double freq) {
  PortMode mode{};
  mode.fc = c0 / (2.0 * port.a);
  const double k = 2.0 * M_PI * freq / c0;
  const double kc = M_PI / port.a; // TE10 cutoff wavenumber
  const double beta = std::sqrt(std::max(0.0, k * k - kc * kc));
  if (beta > 0) {
    mode.Z0 = eta0 * k / beta;
    const double area = port.a * port.b;
    mode.E0 = std::sqrt(4.0 * mode.Z0 / area); // 1 W normalization
  } else {
    mode.Z0 = 0.0;
    mode.E0 = 0.0;
  }
  return mode;
}

SParams2 straight_waveguide_sparams(const RectWaveguidePort &port,
                                    double length, double freq) {
  const double k = 2.0 * M_PI * freq / c0;
  const double kc = M_PI / port.a;
  const double beta = std::sqrt(std::max(0.0, k * k - kc * kc));
  const std::complex<double> phase =
      std::exp(std::complex<double>(0.0, -beta * length));
  SParams2 s{};
  s.s11 = 0.0;
  s.s22 = 0.0;
  s.s21 = phase;
  s.s12 = phase;
  return s;
}

} // namespace vectorem

