#include <cassert>
#include <cstdio>
#include <vector>

#include "vectorem/mor/vector_fit.hpp"
#include "vectorem/sweep.hpp"

using namespace vectorem;
using namespace vectorem::mor;

int main() {
  std::vector<double> freqs;
  const int N = 401;
  for (int i = 0; i < N; ++i)
    freqs.push_back(1.0 + (2.0 - 1.0) * i / (N - 1));

  auto log = sweep_two_resonator(freqs, SweepPolicy::Discrete);
  const char *path = "vector_fit_model.json";
  auto fit = vector_fit(freqs, log.s21, 2, path);
  assert(fit.rms_error < 0.01);
  assert(fit.passive);
  std::FILE *fp = std::fopen(path, "r");
  assert(fp != nullptr);
  std::fclose(fp);
  std::remove(path);
  return 0;
}
