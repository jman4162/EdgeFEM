#pragma once

#include <vector>

#include <Eigen/Core>

#include "edgefem/bc.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"

namespace edgefem {

/// Far-field point in 2D pattern cut.
struct FFPoint2D {
  double theta;      // Angle in radians
  double magnitude;  // Field magnitude (linear, not dB)
  double phase;      // Phase in radians
};

/// Embedded element pattern for one port.
struct EmbeddedPattern {
  int port_index;                       // Index of the active port
  std::vector<FFPoint2D> pattern_phi0;  // E-plane cut (phi=0)
  std::vector<FFPoint2D> pattern_phi90; // H-plane cut (phi=90)
  double frequency;                     // Operating frequency (Hz)
  Eigen::Vector3d element_position;     // Physical position of element
};

/// Configuration for embedded pattern computation.
struct ArrayPatternConfig {
  std::vector<double> theta_rad;       // Theta angles for pattern cuts
  double phi_e_plane = 0.0;            // Phi angle for E-plane (usually 0)
  double phi_h_plane = M_PI / 2.0;     // Phi angle for H-plane (usually π/2)
  double k0;                           // Free-space wavenumber (rad/m)
  double huygens_radius = 0.0;         // Radius for Huygens surface (0 = auto)
};

/// Compute embedded element patterns for all ports.
/// The embedded pattern includes mutual coupling effects from other elements.
///
/// @param mesh Volume mesh
/// @param params Maxwell parameters
/// @param bc Boundary conditions
/// @param ports Wave port definitions (one per element)
/// @param config Pattern configuration
/// @return Vector of EmbeddedPattern, one per port
std::vector<EmbeddedPattern> compute_embedded_patterns(
    const Mesh &mesh, const MaxwellParams &params, const BC &bc,
    const std::vector<WavePort> &ports, const ArrayPatternConfig &config);

/// Compute embedded pattern for a single port.
/// @param mesh Volume mesh
/// @param params Maxwell parameters
/// @param bc Boundary conditions
/// @param ports All ports (for proper loading)
/// @param active_port Index of active port
/// @param config Pattern configuration
/// @return EmbeddedPattern for the active port
EmbeddedPattern compute_single_embedded_pattern(
    const Mesh &mesh, const MaxwellParams &params, const BC &bc,
    const std::vector<WavePort> &ports, int active_port,
    const ArrayPatternConfig &config);

/// Convert linear magnitude to dB.
inline double to_db(double linear) {
  return 20.0 * std::log10(std::max(linear, 1e-20));
}

/// Normalize pattern to peak value of 1.0.
void normalize_pattern(std::vector<FFPoint2D> &pattern);

} // namespace edgefem
