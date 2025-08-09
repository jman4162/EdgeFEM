#pragma once

#include <complex>

namespace vectorem {

struct RectWaveguidePort {
  double a; // width (m)
  double b; // height (m)
};

struct PortMode {
  double fc; // cutoff frequency (Hz)
  double Z0; // modal impedance (ohms)
  double E0; // electric field peak for 1 W (V/m)
};

struct SParams2 {
  std::complex<double> s11;
  std::complex<double> s21;
  std::complex<double> s12;
  std::complex<double> s22;
};

PortMode solve_te10_mode(const RectWaveguidePort &port, double freq);

SParams2 straight_waveguide_sparams(const RectWaveguidePort &port,
                                    double length, double freq);

} // namespace vectorem

