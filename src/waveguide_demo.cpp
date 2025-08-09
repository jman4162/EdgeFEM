#include <cmath>
#include <cstdlib>
#include <iostream>

#include "vectorem/ports/port_eigensolve.hpp"

using namespace vectorem;

int main(int argc, char **argv) {
  double f0 = 5e9;  // Hz
  double f1 = 10e9; // Hz
  int N = 51;
  if (argc == 4) {
    f0 = std::atof(argv[1]) * 1e9;
    f1 = std::atof(argv[2]) * 1e9;
    N = std::atoi(argv[3]);
  }
  const RectWaveguidePort port{0.02286, 0.01016}; // WR-90
  const double length = 0.05;                     // 5 cm section
  std::cout << "freq_hz,s11_mag,s21_mag\n";
  for (int i = 0; i < N; ++i) {
    double f = f0 + (f1 - f0) * i / (N - 1);
    auto s = straight_waveguide_sparams(port, length, f);
    std::cout << f << ',' << std::abs(s.s11) << ',' << std::abs(s.s21) << '\n';
  }
  return 0;
}
