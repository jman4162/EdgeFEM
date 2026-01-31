#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include "vectorem/coupling.hpp"

using namespace vectorem;

void test_compute_coupling_basic() {
  std::cout << "test_compute_coupling_basic..." << std::endl;

  // Simple matched network: S = 0
  Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(2, 2);

  auto cm = compute_coupling(S, 10e9, 50.0);

  assert(cm.num_ports == 2);
  assert(cm.z0 == 50.0);
  assert(cm.frequency == 10e9);

  // For S = 0: Z = z0 * I, Y = I/z0
  Eigen::MatrixXcd expected_Z = 50.0 * Eigen::MatrixXcd::Identity(2, 2);
  assert(cm.Z.isApprox(expected_Z, 1e-10));

  std::cout << "  Z matrix:\n" << cm.Z << std::endl;
}

void test_compute_coupling_short() {
  std::cout << "test_compute_coupling_short..." << std::endl;

  // Short circuit: S11 = -1
  Eigen::MatrixXcd S(1, 1);
  S(0, 0) = std::complex<double>(-1.0, 0.0);

  auto cm = compute_coupling(S, 10e9, 50.0);

  // For S = -1: Z = 0 (short circuit)
  assert(std::abs(cm.Z(0, 0)) < 1e-10);

  std::cout << "  Short circuit Z = " << cm.Z(0, 0) << std::endl;
}

void test_mutual_coupling_db() {
  std::cout << "test_mutual_coupling_db..." << std::endl;

  Eigen::MatrixXcd S(3, 3);
  S << std::complex<double>(0.1, 0.0), std::complex<double>(0.01, 0.0),
      std::complex<double>(0.001, 0.0), std::complex<double>(0.01, 0.0),
      std::complex<double>(0.1, 0.0), std::complex<double>(0.01, 0.0),
      std::complex<double>(0.001, 0.0), std::complex<double>(0.01, 0.0),
      std::complex<double>(0.1, 0.0);

  auto cm = compute_coupling(S, 10e9);
  auto db = mutual_coupling_db(cm);

  // Check diagonal (return loss)
  double s11_db = 20.0 * std::log10(0.1);
  assert(std::abs(db(0, 0) - s11_db) < 0.1);

  // Check off-diagonal
  double s12_db = 20.0 * std::log10(0.01);
  assert(std::abs(db(0, 1) - s12_db) < 0.1);

  std::cout << "  S11 = " << db(0, 0) << " dB (expected " << s11_db << ")"
            << std::endl;
  std::cout << "  S12 = " << db(0, 1) << " dB (expected " << s12_db << ")"
            << std::endl;
}

void test_active_reflection() {
  std::cout << "test_active_reflection..." << std::endl;

  // 2-port network with known coupling
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.1, 0.0), std::complex<double>(0.3, 0.0),
      std::complex<double>(0.3, 0.0), std::complex<double>(0.1, 0.0);

  // Uniform excitation (both ports excited with same amplitude)
  Eigen::VectorXcd excitation(2);
  excitation << std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0);

  auto gamma = active_reflection(S, excitation);

  // Active reflection = (S11*a1 + S12*a2) / a1 = S11 + S12 = 0.4
  std::complex<double> expected(0.4, 0.0);
  assert(std::abs(gamma(0) - expected) < 1e-10);

  std::cout << "  Active Gamma = " << gamma(0) << " (expected " << expected
            << ")" << std::endl;
}

void test_active_reflection_phased() {
  std::cout << "test_active_reflection_phased..." << std::endl;

  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.1, 0.0), std::complex<double>(0.3, 0.0),
      std::complex<double>(0.3, 0.0), std::complex<double>(0.1, 0.0);

  // 180° phase difference between ports
  Eigen::VectorXcd excitation(2);
  excitation << std::complex<double>(1.0, 0.0), std::complex<double>(-1.0, 0.0);

  auto gamma = active_reflection(S, excitation);

  // Active reflection = S11 - S12 = 0.1 - 0.3 = -0.2
  std::complex<double> expected(-0.2, 0.0);
  assert(std::abs(gamma(0) - expected) < 1e-10);

  std::cout << "  180° phased Gamma = " << gamma(0) << std::endl;
}

void test_is_passive_true() {
  std::cout << "test_is_passive_true..." << std::endl;

  // Lossless reciprocal network
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.0, 0.3), std::complex<double>(0.0, 0.9),
      std::complex<double>(0.0, 0.9), std::complex<double>(0.0, 0.3);

  // Check |S11|^2 + |S21|^2 = 0.09 + 0.81 = 0.9 < 1
  assert(is_passive(S));

  std::cout << "  Network is passive: true" << std::endl;
}

void test_is_passive_false() {
  std::cout << "test_is_passive_false..." << std::endl;

  // Active network (gain)
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.0, 0.0), std::complex<double>(2.0, 0.0),
      std::complex<double>(2.0, 0.0), std::complex<double>(0.0, 0.0);

  assert(!is_passive(S));

  std::cout << "  Active network detected correctly" << std::endl;
}

void test_is_lossless() {
  std::cout << "test_is_lossless..." << std::endl;

  // Ideal 3dB splitter (lossless)
  double a = 1.0 / std::sqrt(2.0);
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, a),
      std::complex<double>(0.0, a), std::complex<double>(0.0, 0.0);

  // This is NOT lossless: |S12|^2 + |S22|^2 = 0.5 ≠ 1
  assert(!is_lossless(S, 1e-6));

  // Create truly lossless network
  Eigen::MatrixXcd S_lossless(2, 2);
  S_lossless << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 1.0),
      std::complex<double>(0.0, 1.0), std::complex<double>(0.0, 0.0);

  assert(is_lossless(S_lossless, 1e-6));

  std::cout << "  Lossless check working correctly" << std::endl;
}

void test_is_reciprocal() {
  std::cout << "test_is_reciprocal..." << std::endl;

  // Reciprocal network
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.1, 0.1), std::complex<double>(0.9, 0.0),
      std::complex<double>(0.9, 0.0), std::complex<double>(0.1, 0.1);

  assert(is_reciprocal(S));

  // Non-reciprocal (isolator)
  Eigen::MatrixXcd S_nr(2, 2);
  S_nr << std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0);

  assert(!is_reciprocal(S_nr, 1e-6));

  std::cout << "  Reciprocity check working correctly" << std::endl;
}

int main() {
  test_compute_coupling_basic();
  test_compute_coupling_short();
  test_mutual_coupling_db();
  test_active_reflection();
  test_active_reflection_phased();
  test_is_passive_true();
  test_is_passive_false();
  test_is_lossless();
  test_is_reciprocal();
  std::cout << "All coupling tests passed!" << std::endl;
  return 0;
}
