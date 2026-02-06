#pragma once

#include <complex>
#include <vector>

namespace edgefem {

enum class SweepPolicy { Discrete, Reuse, Balanced };

struct SweepLog {
  std::vector<double> solve_time;        // wall-clock seconds per solve
  std::vector<int> iterations;           // iteration count per solve
  std::vector<std::complex<double>> s21; // second port response
};

// Runs a sweep of a toy two-resonator system over the given frequencies.
// Returns timing and iteration logs as well as the S21 response at each point.
SweepLog sweep_two_resonator(const std::vector<double> &freqs,
                             SweepPolicy policy);

} // namespace edgefem
