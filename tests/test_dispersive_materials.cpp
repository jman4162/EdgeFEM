#include "edgefem/materials/dispersive.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

using namespace edgefem::materials;

constexpr double PI = 3.14159265358979323846;

// Helper to check if two complex numbers are approximately equal
bool approx_equal(std::complex<double> a, std::complex<double> b, double tol = 1e-6) {
  return std::abs(a - b) < tol;
}

bool approx_equal_rel(std::complex<double> a, std::complex<double> b, double rel_tol = 1e-4) {
  double mag = std::max(std::abs(a), std::abs(b));
  if (mag < 1e-15) return true;
  return std::abs(a - b) / mag < rel_tol;
}

void test_debye_basic() {
  std::cout << "Testing Debye model basic functionality..." << std::endl;

  // Simple Debye material
  double eps_s = 80.0;   // static permittivity (like water)
  double eps_inf = 4.0;  // optical permittivity
  double tau = 1e-11;    // relaxation time ~10 ps

  DebyeMaterial debye(eps_s, eps_inf, tau);

  // At DC (omega=0), should return eps_s
  auto eps_dc = debye.eval_eps(0.0);
  assert(approx_equal(eps_dc, std::complex<double>(eps_s, 0.0), 1e-10));
  std::cout << "  DC limit: eps = " << eps_dc << " (expected " << eps_s << ")" << std::endl;

  // At very high frequency, should approach eps_inf
  double omega_high = 1e15;  // ~159 THz
  auto eps_high = debye.eval_eps(omega_high);
  assert(std::abs(eps_high.real() - eps_inf) < 0.1);
  std::cout << "  High-freq limit: eps = " << eps_high << " (expected ~" << eps_inf << ")" << std::endl;

  // At omega*tau = 1, imaginary part should be maximum
  double omega_peak = 1.0 / tau;
  auto eps_peak = debye.eval_eps(omega_peak);
  std::cout << "  At omega*tau=1: eps = " << eps_peak << std::endl;

  // Verify accessor methods
  assert(debye.eps_static() == eps_s);
  assert(debye.eps_inf() == eps_inf);
  assert(debye.tau() == tau);

  std::cout << "  PASSED" << std::endl;
}

void test_debye_water() {
  std::cout << "Testing Debye model for water at 20°C..." << std::endl;

  // Water at 20°C: eps_s = 80.1, eps_inf = 4.9, tau = 9.3 ps
  DebyeMaterial water(80.1, 4.9, 9.3e-12);

  // Test at 10 GHz
  double f_ghz = 10.0;
  double omega = 2.0 * PI * f_ghz * 1e9;
  auto eps = water.eval_eps(omega);

  // Expected values from literature (approximate)
  // At 10 GHz: eps' ~ 55, eps'' ~ 35
  std::cout << "  At 10 GHz: eps = " << eps.real() << " - j" << -eps.imag() << std::endl;
  assert(eps.real() > 40 && eps.real() < 70);  // Reasonable range
  assert(eps.imag() < 0);  // Loss should be negative imaginary (eps'' > 0)

  std::cout << "  PASSED" << std::endl;
}

void test_lorentz_single_pole() {
  std::cout << "Testing Lorentz model single pole..." << std::endl;

  LorentzMaterial lorentz(1.0);  // eps_inf = 1
  double omega0 = 2.0 * PI * 10e9;  // 10 GHz resonance
  double gamma = 2.0 * PI * 0.5e9;  // 0.5 GHz damping
  double delta_eps = 5.0;

  lorentz.add_pole(delta_eps, omega0, gamma);

  // Far below resonance: eps should approach eps_inf + delta_eps
  double omega_low = omega0 / 100.0;
  auto eps_low = lorentz.eval_eps(omega_low);
  std::cout << "  Below resonance: eps = " << eps_low << std::endl;
  assert(eps_low.real() > delta_eps);  // Should be above delta_eps

  // At resonance: maximum imaginary part
  auto eps_res = lorentz.eval_eps(omega0);
  std::cout << "  At resonance: eps = " << eps_res << std::endl;
  assert(std::abs(eps_res.imag()) > 0.1);  // Should have significant loss

  // Far above resonance: eps should approach eps_inf
  double omega_high = omega0 * 100.0;
  auto eps_high = lorentz.eval_eps(omega_high);
  std::cout << "  Above resonance: eps = " << eps_high << std::endl;
  assert(std::abs(eps_high.real() - 1.0) < 0.1);  // Should approach eps_inf

  // Verify structure
  assert(lorentz.num_poles() == 1);
  assert(lorentz.eps_inf() == 1.0);

  std::cout << "  PASSED" << std::endl;
}

void test_lorentz_multi_pole() {
  std::cout << "Testing Lorentz model multi-pole..." << std::endl;

  LorentzMaterial lorentz(2.0);  // eps_inf = 2

  // Add two poles
  lorentz.add_pole(3.0, 2.0 * PI * 5e9, 2.0 * PI * 0.2e9);   // 5 GHz
  lorentz.add_pole(2.0, 2.0 * PI * 15e9, 2.0 * PI * 0.3e9);  // 15 GHz

  assert(lorentz.num_poles() == 2);

  // Test at various frequencies
  auto eps_1ghz = lorentz.eval_eps(2.0 * PI * 1e9);
  auto eps_10ghz = lorentz.eval_eps(2.0 * PI * 10e9);
  auto eps_50ghz = lorentz.eval_eps(2.0 * PI * 50e9);

  std::cout << "  At 1 GHz: eps = " << eps_1ghz << std::endl;
  std::cout << "  At 10 GHz: eps = " << eps_10ghz << std::endl;
  std::cout << "  At 50 GHz: eps = " << eps_50ghz << std::endl;

  // Low frequency should have high eps (all poles contribute)
  assert(eps_1ghz.real() > 5.0);
  // High frequency should approach eps_inf
  assert(std::abs(eps_50ghz.real() - 2.0) < 0.5);

  std::cout << "  PASSED" << std::endl;
}

void test_drude_basic() {
  std::cout << "Testing Drude model basic functionality..." << std::endl;

  double omega_p = 1.37e16;  // Plasma frequency (gold-like)
  double gamma = 4.05e13;    // Damping (gold-like)

  DrudeMaterial drude(omega_p, gamma);

  // Below plasma frequency: eps should be negative (metallic)
  double omega_low = omega_p / 10.0;
  auto eps_low = drude.eval_eps(omega_low);
  std::cout << "  At omega_p/10: eps = " << eps_low << std::endl;
  assert(eps_low.real() < 0);  // Metallic behavior

  // At plasma frequency: eps_real should be near zero
  auto eps_plasma = drude.eval_eps(omega_p);
  std::cout << "  At omega_p: eps = " << eps_plasma << std::endl;
  // Real part should be small compared to 1
  assert(std::abs(eps_plasma.real()) < 1.0);

  // Above plasma frequency: eps should approach 1
  double omega_high = omega_p * 10.0;
  auto eps_high = drude.eval_eps(omega_high);
  std::cout << "  At 10*omega_p: eps = " << eps_high << std::endl;
  assert(eps_high.real() > 0.9 && eps_high.real() < 1.1);

  // Verify accessor methods
  assert(drude.omega_p() == omega_p);
  assert(drude.gamma() == gamma);

  std::cout << "  PASSED" << std::endl;
}

void test_drude_gold_optical() {
  std::cout << "Testing Drude model for gold at optical frequencies..." << std::endl;

  // Gold parameters (approximate)
  double omega_p = 1.37e16;  // rad/s
  double gamma = 4.05e13;    // rad/s

  DrudeMaterial gold(omega_p, gamma);

  // Test at 500 nm wavelength (visible light)
  double c = 2.998e8;
  double lambda = 500e-9;
  double omega = 2.0 * PI * c / lambda;

  auto eps = gold.eval_eps(omega);
  std::cout << "  At 500 nm: eps = " << eps.real() << " + j" << eps.imag() << std::endl;

  // Gold at visible frequencies has negative real part
  assert(eps.real() < 0);

  // Test at 1550 nm (telecom wavelength)
  lambda = 1550e-9;
  omega = 2.0 * PI * c / lambda;
  eps = gold.eval_eps(omega);
  std::cout << "  At 1550 nm: eps = " << eps.real() << " + j" << eps.imag() << std::endl;
  assert(eps.real() < 0);  // Still metallic

  std::cout << "  PASSED" << std::endl;
}

void test_drude_lorentz_combined() {
  std::cout << "Testing Drude-Lorentz combined model..." << std::endl;

  // Gold with interband transitions
  double eps_inf = 1.0;
  double omega_p = 1.37e16;
  double gamma_d = 4.05e13;

  DrudeLorentzMaterial gold_dl(eps_inf, omega_p, gamma_d);

  // Add interband transition (simplified)
  double omega_lorentz = 4.0e15;  // ~0.5 eV
  gold_dl.add_lorentz_pole(2.0, omega_lorentz, 1e14);

  // Test at optical frequency
  double c = 2.998e8;
  double omega = 2.0 * PI * c / 600e-9;  // 600 nm

  auto eps = gold_dl.eval_eps(omega);
  std::cout << "  At 600 nm: eps = " << eps.real() << " + j" << eps.imag() << std::endl;

  // Should still be metallic but modified by Lorentz pole
  assert(eps.real() < 5);  // Negative or small positive due to interband

  // Verify structure
  assert(gold_dl.lorentz_poles().size() == 1);
  assert(gold_dl.omega_p() == omega_p);

  std::cout << "  PASSED" << std::endl;
}

void test_parameter_validation() {
  std::cout << "Testing parameter validation..." << std::endl;

  // Debye: tau must be positive
  bool caught = false;
  try {
    DebyeMaterial bad_debye(80.0, 4.0, -1e-12);
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  assert(caught);
  std::cout << "  Debye negative tau: correctly rejected" << std::endl;

  // Drude: omega_p must be positive
  caught = false;
  try {
    DrudeMaterial bad_drude(0.0, 1e13);
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  assert(caught);
  std::cout << "  Drude zero omega_p: correctly rejected" << std::endl;

  // Drude: gamma must be non-negative
  caught = false;
  try {
    DrudeMaterial bad_drude2(1e16, -1e13);
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  assert(caught);
  std::cout << "  Drude negative gamma: correctly rejected" << std::endl;

  // Lorentz: omega0 must be positive
  caught = false;
  try {
    LorentzMaterial lorentz;
    lorentz.add_pole(1.0, -1e9, 1e8);
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  assert(caught);
  std::cout << "  Lorentz negative omega0: correctly rejected" << std::endl;

  // Lorentz: gamma must be non-negative
  caught = false;
  try {
    LorentzMaterial lorentz;
    lorentz.add_pole(1.0, 1e9, -1e8);
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  assert(caught);
  std::cout << "  Lorentz negative gamma: correctly rejected" << std::endl;

  std::cout << "  PASSED" << std::endl;
}

void test_mu_default() {
  std::cout << "Testing default mu_r behavior..." << std::endl;

  DebyeMaterial debye(80.0, 4.0, 1e-11);
  DrudeMaterial drude(1e16, 1e13);
  LorentzMaterial lorentz;
  lorentz.add_pole(1.0, 1e10, 1e8);

  // All non-magnetic materials should return mu_r = 1
  double omega = 2.0 * PI * 10e9;

  assert(approx_equal(debye.eval_mu(omega), std::complex<double>(1.0, 0.0)));
  assert(approx_equal(drude.eval_mu(omega), std::complex<double>(1.0, 0.0)));
  assert(approx_equal(lorentz.eval_mu(omega), std::complex<double>(1.0, 0.0)));

  std::cout << "  All materials return mu_r = 1.0" << std::endl;
  std::cout << "  PASSED" << std::endl;
}

void test_silicon_10ghz() {
  std::cout << "Testing silicon at 10 GHz (verification case)..." << std::endl;

  // Silicon: eps_s ~ 11.68, eps_inf ~ 11.0, tau ~ 1e-13 s
  // At 10 GHz, the dispersion is minimal (omega*tau << 1)
  DebyeMaterial silicon(11.68, 11.0, 1e-13);

  double freq = 10e9;
  double omega = 2.0 * PI * freq;
  auto eps = silicon.eval_eps(omega);

  std::cout << "  At 10 GHz: eps = " << eps.real() << " - j" << -eps.imag() << std::endl;

  // At 10 GHz with tau=1e-13, omega*tau ~ 6e-3 << 1
  // So eps should be very close to eps_s
  assert(std::abs(eps.real() - 11.68) < 0.1);
  assert(std::abs(eps.imag()) < 0.1);

  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "=== Dispersive Materials Unit Tests ===" << std::endl;
  std::cout << std::endl;

  test_debye_basic();
  test_debye_water();
  test_lorentz_single_pole();
  test_lorentz_multi_pole();
  test_drude_basic();
  test_drude_gold_optical();
  test_drude_lorentz_combined();
  test_parameter_validation();
  test_mu_default();
  test_silicon_10ghz();

  std::cout << std::endl;
  std::cout << "=== All tests PASSED ===" << std::endl;
  return 0;
}
