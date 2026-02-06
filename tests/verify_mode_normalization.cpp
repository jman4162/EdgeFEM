// Verify TE10 mode field normalization against analytical formula
// The mode field should be normalized for unit power flow.

#include <cmath>
#include <iostream>
#include <iomanip>

int main() {
  std::cout << std::setprecision(10);

  // Physical constants
  constexpr double c0 = 299792458.0;
  constexpr double mu0 = 4.0 * M_PI * 1e-7;
  constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

  // WR-90 at 10 GHz
  double a = 0.02286;  // width
  double b = 0.01016;  // height
  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k = omega / c0;
  double kc = M_PI / a;  // TE10 cutoff wavenumber
  double beta = std::sqrt(k*k - kc*kc);

  std::cout << "=== TE10 Mode Normalization Check ===" << std::endl;
  std::cout << "WR-90: a = " << a*1000 << " mm, b = " << b*1000 << " mm" << std::endl;
  std::cout << "freq = " << freq/1e9 << " GHz" << std::endl;
  std::cout << "omega = " << omega << " rad/s" << std::endl;
  std::cout << "kc = " << kc << " rad/m" << std::endl;
  std::cout << "beta = " << beta << " rad/m" << std::endl;
  std::cout << std::endl;

  // TE10 mode fields (Pozar convention):
  // Hz = A cos(πx/a)
  // Ey = (jωμ₀/kc²)(π/a) A sin(πx/a)
  // Hx = (jβ/kc²)(π/a) A sin(πx/a)
  //
  // Time-averaged Poynting vector:
  // <Sz> = ½ Re{Ey Hx*}
  //      = ½ Re{(jωμ₀π/akc²)(jβπ/akc²)* A² sin²}
  //      = ½ Re{(jωμ₀π/akc²)(-jβπ/akc²) A² sin²}
  //      = ½ ωμ₀βπ²/(a²kc⁴) A² sin²
  //
  // Total power:
  // P = ∫₀ᵃ∫₀ᵇ <Sz> dxdy
  //   = ½ ωμ₀βπ²/(a²kc⁴) A² ∫sin²(πx/a)dx ∫dy
  //   = ½ ωμ₀βπ²/(a²kc⁴) A² (a/2) b
  //   = ωμ₀βπ²bA² / (4a kc⁴)
  //
  // Using kc = π/a, so kc⁴ = π⁴/a⁴:
  // P = ωμ₀βπ²bA² / (4a × π⁴/a⁴)
  //   = ωμ₀βπ²bA² × a⁴ / (4a × π⁴)
  //   = ωμ₀βbA²a³ / (4π²)
  //
  // For unit power (P = 1):
  // A² = 4π² / (ωμ₀βba³)

  double pi_sq = M_PI * M_PI;
  double A_sq_correct = 4.0 * pi_sq / (omega * mu0 * beta * b * a*a*a);
  double A_correct = std::sqrt(A_sq_correct);

  std::cout << "=== Correct Normalization ===" << std::endl;
  std::cout << "A² = 4π² / (ωμ₀βba³)" << std::endl;
  std::cout << "A² = " << A_sq_correct << std::endl;
  std::cout << "A = " << A_correct << std::endl;
  std::cout << std::endl;

  // What the code currently uses:
  // A² = 8 * kc⁴ / (ω * μ * β * π² * b)
  double A_sq_code = 8.0 * std::pow(kc, 4) / (omega * mu0 * beta * pi_sq * b);
  double A_code = std::sqrt(A_sq_code);

  std::cout << "=== Current Code ===" << std::endl;
  std::cout << "A² = 8 * kc⁴ / (ωμβπ²b)" << std::endl;
  std::cout << "A² = " << A_sq_code << std::endl;
  std::cout << "A = " << A_code << std::endl;
  std::cout << std::endl;

  // Ratio
  double ratio = A_code / A_correct;
  std::cout << "=== Comparison ===" << std::endl;
  std::cout << "A_code / A_correct = " << ratio << std::endl;
  std::cout << "A²_code / A²_correct = " << ratio*ratio << std::endl;
  std::cout << std::endl;

  // The issue: code formula is missing factor of 'a' in numerator
  // Code: A² = 8 * kc⁴ / (ωμβπ²b)
  // With kc⁴ = π⁴/a⁴:
  //       A² = 8π⁴ / (a⁴ωμβπ²b) = 8π² / (a⁴ωμβb)
  //
  // Correct: A² = 4π² / (a³ωμβb)
  //
  // Ratio: (8π²/a⁴) / (4π²/a³) = 2/a
  std::cout << "Expected ratio² = 2/a = " << 2.0/a << std::endl;
  std::cout << "Actual ratio² = " << ratio*ratio << std::endl;
  std::cout << std::endl;

  // Verify power integral with correct normalization
  double power_factor = omega * mu0 * beta * b / (4.0 * pi_sq);  // × A² × a³ gives P
  double P_correct = power_factor * A_sq_correct * a*a*a;
  double P_code = power_factor * A_sq_code * a*a*a;

  std::cout << "=== Power Verification ===" << std::endl;
  std::cout << "Power with A_correct: " << P_correct << " (should be 1)" << std::endl;
  std::cout << "Power with A_code: " << P_code << " (should be 1)" << std::endl;
  std::cout << std::endl;

  // The fix: multiply by 'a' or equivalently divide by kc/π
  // New formula: A² = 4π² / (a³ωμβb)
  // Or: A² = 4 * kc⁴ * a / (ωμβπ²b)
  double A_sq_fixed = 4.0 * std::pow(kc, 4) * a / (omega * mu0 * beta * pi_sq * b);
  double A_fixed = std::sqrt(A_sq_fixed);
  double P_fixed = power_factor * A_sq_fixed * a*a*a;

  std::cout << "=== Fixed Formula ===" << std::endl;
  std::cout << "A² = 4 * kc⁴ * a / (ωμβπ²b)" << std::endl;
  std::cout << "A² = " << A_sq_fixed << std::endl;
  std::cout << "A = " << A_fixed << std::endl;
  std::cout << "Power with A_fixed: " << P_fixed << " (should be 1)" << std::endl;
  std::cout << std::endl;

  // Check ratio
  std::cout << "=== Final Check ===" << std::endl;
  std::cout << "A_fixed / A_correct = " << A_fixed / A_correct << " (should be 1)" << std::endl;

  // Impact on port weights
  // Port weight ~ Et × edge_length ~ (ωμ/kc²) × A × (1/a) × edge_len
  // ||w||² ~ A² × (constant geometric factor)
  // So ||w||²_code / ||w||²_correct = A²_code / A²_correct = 87.5
  std::cout << std::endl;
  std::cout << "=== Impact on Port Weights ===" << std::endl;
  std::cout << "||w||²_code / ||w||²_correct = " << A_sq_code / A_sq_correct << std::endl;
  std::cout << "This means port admittance Y = ||w||²/Z0 is " << A_sq_code / A_sq_correct << "x too large!" << std::endl;

  return 0;
}
