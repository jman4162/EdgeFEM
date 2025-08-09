#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "vectorem/sweep.hpp"

using namespace vectorem;

int main(int argc, char **argv) {
  SweepPolicy policy = SweepPolicy::Discrete;
  if (argc > 1) {
    std::string arg = argv[1];
    std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
    if (arg == "reuse")
      policy = SweepPolicy::Reuse;
    else if (arg == "balanced")
      policy = SweepPolicy::Balanced;
  }

  std::vector<double> freqs;
  const int N = 401;
  freqs.reserve(N);
  for (int i = 0; i < N; ++i)
    freqs.push_back(1.0 + (2.0 - 1.0) * i / (N - 1));

  auto log = sweep_two_resonator(freqs, policy);
  double total =
      std::accumulate(log.solve_time.begin(), log.solve_time.end(), 0.0);
  std::cout << "policy=" << static_cast<int>(policy) << " total_time=" << total
            << " s\n";
  return 0;
}
