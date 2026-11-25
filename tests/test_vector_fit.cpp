#include "test_check.hpp"
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "edgefem/mor/vector_fit.hpp"
#include "edgefem/sweep.hpp"

using namespace edgefem;
using namespace edgefem::mor;

int main() {
  std::vector<double> freqs;
  const int N = 401;
  for (int i = 0; i < N; ++i)
    freqs.push_back(1.0 + (2.0 - 1.0) * i / (N - 1));

  auto log = sweep_two_resonator(freqs, SweepPolicy::Discrete);
  const char *path = "vector_fit_model.json";
  auto fit = vector_fit(freqs, log.s21, 2, path);
  CHECK(fit.rms_error < 0.01);
  CHECK(fit.passive);
  std::FILE *fp = std::fopen(path, "r");
  CHECK(fp != nullptr);
  std::fclose(fp);
  std::remove(path);

  // Inputs the least-squares system cannot determine are refused: fewer
  // samples than the 2n + 1 unknowns, mismatched lengths, order < 1.
  auto throws = [](auto &&f) {
    try {
      f();
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  const std::vector<double> f4{1.0, 1.1, 1.2, 1.3};
  const std::vector<std::complex<double>> s4(4, {0.5, 0.1});
  CHECK(throws([&] { vector_fit(f4, s4, 2, path); })); // 4 < 5 unknowns
  CHECK(throws([&] { vector_fit(f4, {s4[0], s4[1]}, 1, path); }));
  CHECK(throws([&] { vector_fit(f4, s4, 0, path); }));
  std::remove(path);
  return edgefem_test::exit_code();
}
