#pragma once

#include <Eigen/Core>

#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"

namespace edgefem {

struct LumpedPortConfig {
  int surface_tag;                            // Physical tag of port gap surface
  double z0 = 50.0;                           // Reference impedance (ohms)
  Eigen::Vector3d e_direction = {0, 0, 1};    // E-field direction across gap
};

/// Build a lumped port from a tagged surface.
/// The returned WavePort has weights computed by projecting edge vectors onto
/// the specified E-field direction. It plugs directly into assemble_maxwell()
/// with no modifications — the assembly only uses edges, weights, and mode.Z0.
///
/// @param mesh Volume mesh
/// @param config Lumped port configuration (surface tag, Z0, E-field direction)
/// @return WavePort compatible with assemble_maxwell()
WavePort build_lumped_port(const Mesh &mesh, const LumpedPortConfig &config);

} // namespace edgefem
