#include "edgefem/mor/vector_fit.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <sstream>

namespace edgefem {
namespace mor {

FitResult vector_fit(const std::vector<double> &freq,
                     const std::vector<std::complex<double>> &sparam, int order,
                     const std::string &out_path) {
  const size_t N = freq.size();
  const int cols = 2 * order + 1; // a1..an, b0..bn
  Eigen::MatrixXcd A(N, cols);
  Eigen::VectorXcd y(N);

  for (size_t i = 0; i < N; ++i) {
    std::complex<double> s = sparam[i];
    double f = freq[i];
    y[i] = s;
    for (int k = 0; k < order; ++k) {
      A(i, k) = -s * std::pow(f, k + 1);
    }
    A(i, order) = 1.0;
    for (int k = 0; k < order; ++k) {
      A(i, order + 1 + k) = std::pow(f, k + 1);
    }
  }

  Eigen::VectorXcd sol = A.colPivHouseholderQr().solve(y);
  FitResult res;
  res.denominator.push_back({1.0, 0.0});
  for (int k = 0; k < order; ++k)
    res.denominator.push_back(sol[k]);
  res.numerator.push_back(sol[order]);
  for (int k = 0; k < order; ++k)
    res.numerator.push_back(sol[order + 1 + k]);

  std::vector<std::complex<double>> fit(N);
  for (size_t i = 0; i < N; ++i) {
    std::complex<double> num = res.numerator[0];
    std::complex<double> den = res.denominator[0];
    double f = freq[i];
    for (int k = 1; k < res.numerator.size(); ++k)
      num += res.numerator[k] * std::pow(f, k);
    for (int k = 1; k < res.denominator.size(); ++k)
      den += res.denominator[k] * std::pow(f, k);
    fit[i] = num / den;
  }

  double err = 0.0;
  double ref = 0.0;
  double max_mag = 0.0;
  for (size_t i = 0; i < N; ++i) {
    err += std::norm(fit[i] - sparam[i]);
    ref += std::norm(sparam[i]);
    max_mag = std::max(max_mag, std::abs(fit[i]));
  }
  res.rms_error = std::sqrt(err / N) / std::sqrt(ref / N);
  res.passive = max_mag <= 1.0 + 1e-6;

  std::ofstream os(out_path);
  os << "{\n";
  os << "  \"numerator\": [";
  for (size_t i = 0; i < res.numerator.size(); ++i) {
    os << "[" << res.numerator[i].real() << "," << res.numerator[i].imag()
       << "]";
    if (i + 1 != res.numerator.size())
      os << ",";
  }
  os << "],\n";
  os << "  \"denominator\": [";
  for (size_t i = 0; i < res.denominator.size(); ++i) {
    os << "[" << res.denominator[i].real() << "," << res.denominator[i].imag()
       << "]";
    if (i + 1 != res.denominator.size())
      os << ",";
  }
  os << "],\n";
  os << "  \"rms_error\": " << res.rms_error << ",\n";
  os << "  \"passive\": " << (res.passive ? "true" : "false") << "\n";
  os << "}\n";
  return res;
}

} // namespace mor
} // namespace edgefem
