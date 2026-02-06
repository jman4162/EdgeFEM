#include "edgefem/sweep.hpp"

#include <Eigen/Dense>
#include <chrono>

namespace edgefem {
namespace {

Eigen::Matrix2d build_A(double f) {
  const double f0 = 1.0; // GHz normalized
  Eigen::Matrix2d A0;
  A0 << 4.0, 1.0, 1.0, 3.0;
  Eigen::Matrix2d A1;
  A1 << 0.1, 0.05, 0.05, 0.2;
  return A0 + (f - f0) * A1;
}

// Simple Richardson iteration solver for Ax=b
Eigen::Vector2d solve_iterative(const Eigen::Matrix2d &A,
                                const Eigen::Vector2d &b,
                                const Eigen::Vector2d &x0, int &iters) {
  Eigen::Vector2d x = x0;
  const double omega = 0.1;
  iters = 0;
  for (; iters < 1000; ++iters) {
    Eigen::Vector2d r = b - A * x;
    if (r.norm() < 1e-8)
      break;
    x += omega * r;
  }
  return x;
}

} // namespace

SweepLog sweep_two_resonator(const std::vector<double> &freqs,
                             SweepPolicy policy) {
  SweepLog log;
  log.solve_time.reserve(freqs.size());
  log.iterations.reserve(freqs.size());
  log.s21.reserve(freqs.size());

  Eigen::Vector2d x_prev = Eigen::Vector2d::Zero();
  const Eigen::Vector2d b(1.0, 0.0);
  size_t reset_period = freqs.size() / 10 + 1;

  for (size_t i = 0; i < freqs.size(); ++i) {
    const double f = freqs[i];
    Eigen::Vector2d guess = Eigen::Vector2d::Zero();
    if (policy == SweepPolicy::Reuse && i > 0) {
      guess = x_prev;
    } else if (policy == SweepPolicy::Balanced && i > 0 &&
               (i % reset_period) != 0) {
      guess = x_prev;
    }

    auto A = build_A(f);
    int iters = 0;
    auto start = std::chrono::steady_clock::now();
    auto x = solve_iterative(A, b, guess, iters);
    auto end = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(end - start).count();

    log.solve_time.push_back(dt);
    log.iterations.push_back(iters);
    log.s21.push_back(std::complex<double>(x[1], 0.0));
    x_prev = x;
  }
  return log;
}

} // namespace edgefem
