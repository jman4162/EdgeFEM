#pragma once

#include <vector>

#include <Eigen/Core>

#include "vectorem/mesh.hpp"
#include "vectorem/ports/port_eigensolve.hpp"

namespace vectorem {

struct PortSurfaceMesh {
  Mesh mesh;
  std::vector<int> volume_tri_indices; // index into volume mesh tris
};

PortSurfaceMesh extract_surface_mesh(const Mesh &volume_mesh, int surface_tag);

struct WavePort {
  int surface_tag = 0;
  PortMode mode;
  std::vector<int> edges;      // global edge indices participating in the port
  Eigen::VectorXcd weights;    // line-integral weights for each edge
};

WavePort build_wave_port(const Mesh &volume_mesh, const PortSurfaceMesh &surface,
                         const PortMode &mode);

} // namespace vectorem

