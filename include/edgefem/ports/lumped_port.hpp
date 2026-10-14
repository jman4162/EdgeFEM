#pragma once

#include <Eigen/Core>

#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"

namespace edgefem {

/// Weight computation strategy for lumped ports.
enum class LumpedPortWeightMode {
  /// Projection: w_i = edge_vec . e_dir, scaled so that a uniform gap field
  /// gives the exact voltage. Not mesh-convergent: weight accuracy depends
  /// on element shape.
  Projection,

  /// SurfaceIntegral: w_i = ∫(N_i · e_dir) dS using Whitney basis functions
  /// on the port surface triangles. Mesh-convergent and physically motivated.
  SurfaceIntegral,
};

struct LumpedPortConfig {
  int surface_tag;                         // Physical tag of port gap surface
  double z0 = 50.0;                        // Reference impedance (ohms)
  Eigen::Vector3d e_direction = {0, 0, 1}; // E-field direction across gap
  LumpedPortWeightMode weight_mode = LumpedPortWeightMode::SurfaceIntegral;
};

/// Build a lumped (gap) port from a tagged surface strip.
/// The returned weights are the port voltage functional v = w / W, where
/// w_i = int_S N_i . e_dir dS and W = strip area / (strip extent along
/// e_dir) is the strip width, so the gap voltage is V = v^T x. Use the port
/// directly with assemble_maxwell(), calculate_sparams() and
/// frequency_sweep(); no weight normalization is needed. Throws if e_dir has
/// no extent in the surface (e_dir normal to the strip).
///
/// @param mesh Volume mesh
/// @param config Lumped port configuration (surface tag, Z0, E-field direction)
/// @return WavePort compatible with assemble_maxwell()
WavePort build_lumped_port(const Mesh &mesh, const LumpedPortConfig &config);

} // namespace edgefem
