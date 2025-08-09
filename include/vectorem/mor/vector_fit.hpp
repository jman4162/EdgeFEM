#pragma once

#include <complex>
#include <string>
#include <vector>

namespace vectorem {
namespace mor {

struct FitResult {
  std::vector<std::complex<double>> numerator;   // b0, b1, ...
  std::vector<std::complex<double>> denominator; // 1, a1, ...
  double rms_error = 0.0;
  bool passive = false;
};

// Fit a rational model of given order to S(f) data and write JSON model.
FitResult vector_fit(const std::vector<double> &freq,
                     const std::vector<std::complex<double>> &sparam, int order,
                     const std::string &out_path);

} // namespace mor
} // namespace vectorem
