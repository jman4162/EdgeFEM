#pragma once

#include <complex>
#include <string>
#include <vector>

namespace vectorem {
namespace mor {

/// Coefficients of a fitted rational model \f$H(f)\f$ with a monic denominator.
/// For order *n* the model has the form
/// \f[
///   H(f) = \frac{b_0 + b_1 f + \dots + b_n f^n}{1 + a_1 f + \dots + a_n f^n}
/// \f]
/// where `numerator` holds \f$b_0\ldots b_n\f$ and `denominator` holds
/// \f$1, a_1\ldots a_n\f$.
struct FitResult {
  std::vector<std::complex<double>> numerator;   ///< b0, b1, ...
  std::vector<std::complex<double>> denominator; ///< 1, a1, ...
  double rms_error = 0.0;
  bool passive = false;
};

/// Fit a rational model of a given order to sampled S(f) data and write the
/// coefficients to a JSON file.
///
/// The `order` parameter sets the degree *n* of both the numerator and the
/// denominator as shown above. Because the implementation forms powers of the
/// frequency samples directly, the values in `freq` should be scaled so that
/// \f$|f|\lesssim 1\f$ (e.g. express frequencies in GHz and subtract the
/// mean) to avoid ill-conditioned systems.
FitResult vector_fit(const std::vector<double> &freq,
                     const std::vector<std::complex<double>> &sparam, int order,
                     const std::string &out_path);

} // namespace mor
} // namespace vectorem
