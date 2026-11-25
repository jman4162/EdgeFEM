#include "edgefem/array/embedded_pattern.hpp"

#include <stdexcept>
#include <string>

namespace edgefem {

namespace {
[[noreturn]] void not_implemented(const char *name) {
  throw std::logic_error(
      std::string(name) +
      ": not implemented. The former body summed the raw edge DOFs of every "
      "mesh edge with an angular weight, not a near-to-far transform, so "
      "its patterns had no physical meaning. Compute a pattern with "
      "calculate_sparams_with_fields (or calculate_sparams_eigenmode_with_"
      "fields), extract_huygens_surface and stratton_chu_3d.");
}
} // namespace

std::vector<EmbeddedPattern>
compute_embedded_patterns(const Mesh &, const MaxwellParams &, const BC &,
                          const std::vector<WavePort> &,
                          const ArrayPatternConfig &) {
  not_implemented("compute_embedded_patterns");
}

EmbeddedPattern
compute_single_embedded_pattern(const Mesh &, const MaxwellParams &, const BC &,
                                const std::vector<WavePort> &, int,
                                const ArrayPatternConfig &) {
  not_implemented("compute_single_embedded_pattern");
}

void normalize_pattern(std::vector<FFPoint2D> &pattern) {
  if (pattern.empty())
    return;

  // Find peak magnitude
  double peak = 0.0;
  for (const auto &pt : pattern) {
    if (pt.magnitude > peak) {
      peak = pt.magnitude;
    }
  }

  // Normalize to peak = 1.0
  if (peak > 1e-20) {
    for (auto &pt : pattern) {
      pt.magnitude /= peak;
    }
  }
}

} // namespace edgefem
