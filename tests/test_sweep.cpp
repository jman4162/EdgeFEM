#include "test_check.hpp"
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
  auto bal = sweep_two_resonator(freqs, SweepPolicy::Balanced);
  int it_bal = std::accumulate(bal.iterations.begin(), bal.iterations.end(), 0);
  // Warm starts save iterations on this toy problem: 6464 -> 3563 (Reuse)
  // and 3823 (Balanced, which restarts cold every tenth point).
  CHECK(10 * it_reuse < 6 * it_disc);
  CHECK(it_reuse <= it_bal && it_bal < it_disc);
  return edgefem_test::exit_code();
}
