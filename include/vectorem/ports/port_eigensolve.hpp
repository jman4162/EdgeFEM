#pragma once

#include <complex>
#include <vector>

#include <Eigen/Core>

#include "vectorem/mesh.hpp"

namespace vectorem {

struct RectWaveguidePort {
  double a; // width (m)
  double b; // height (m)
};

enum class ModePolarization { TE, TM };

struct PortMode {
  ModePolarization pol = ModePolarization::TE;
  double fc = 0.0;                         // cutoff frequency (Hz)
  double kc = 0.0;                         // cutoff wavenumber (rad/m)
  double omega = 0.0;                      // operating angular frequency
  std::complex<double> eps = 0.0;          // material permittivity
  std::complex<double> mu = 0.0;           // material permeability
  std::complex<double> beta = 0.0;         // propagation constant
  std::complex<double> Z0 = 0.0;           // modal impedance (ohms)
  Eigen::VectorXcd field;                  // scalar modal field samples
};

struct SParams2 {
  std::complex<double> s11;
  std::complex<double> s21;
  std::complex<double> s12;
  std::complex<double> s22;
};

// Solve for the eigenmodes of a 2D port cross-section.
// Returns a vector of modes, sorted by cutoff frequency.
std::vector<PortMode> solve_port_eigens(
    const Mesh &mesh, int num_modes, double omega,
    std::complex<double> eps_r, std::complex<double> mu_r,
    ModePolarization pol);

PortMode solve_te10_mode(const RectWaveguidePort &port, double freq);

SParams2 straight_waveguide_sparams(const RectWaveguidePort &port,
                                    double length, double freq);

} // namespace vectorem
