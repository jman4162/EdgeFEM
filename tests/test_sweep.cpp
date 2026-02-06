#include <cassert>
#include <numeric>
#include <vector>

#include "edgefem/sweep.hpp"

using namespace edgefem;

int main() {
  std::vector<double> freqs;
  for (int i = 0; i < 101; ++i)
    freqs.push_back(1.0 + (2.0 - 1.0) * i / 100.0);

  auto disc = sweep_two_resonator(freqs, SweepPolicy::Discrete);
  auto reuse = sweep_two_resonator(freqs, SweepPolicy::Reuse);

  int it_disc =
      std::accumulate(disc.iterations.begin(), disc.iterations.end(), 0);
  int it_reuse =
      std::accumulate(reuse.iterations.begin(), reuse.iterations.end(), 0);
  assert(it_disc >= 2 * it_reuse);
  return 0;
}
