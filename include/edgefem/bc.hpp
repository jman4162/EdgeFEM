#pragma once

#include <set>
#include <unordered_set>

#include "edgefem/mesh.hpp"

namespace edgefem {

struct BC {
  std::unordered_set<int> dirichlet_nodes;
  std::unordered_set<int> dirichlet_edges;
};

/// Build scalar PEC boundary condition from mesh surface tag.
/// Returns BC with dirichlet_nodes set for all nodes on surfaces with the given tag.
BC build_scalar_pec(const Mesh &mesh, int pec_tag);

/// Build edge-element PEC boundary condition from mesh surface tag.
/// Returns BC with dirichlet_edges set for all edges on surfaces with the given tag.
/// @note Prints a warning to stderr if no edges are found (common user error).
BC build_edge_pec(const Mesh &mesh, int pec_tag);

/// List all unique physical tags found in the mesh.
/// Returns separate sets for volume (tet) and surface (tri) tags.
/// Useful for debugging when wrong tags are specified.
struct PhysicalTagInfo {
  std::set<int> volume_tags;   ///< Tags from tetrahedra (3D regions)
  std::set<int> surface_tags;  ///< Tags from triangles (2D surfaces)
};
PhysicalTagInfo list_physical_tags(const Mesh &mesh);

/// Check if a surface tag exists in the mesh.
/// @return true if at least one triangle has the given physical tag.
bool has_surface_tag(const Mesh &mesh, int tag);

/// Check if a volume tag exists in the mesh.
/// @return true if at least one tetrahedron has the given physical tag.
bool has_volume_tag(const Mesh &mesh, int tag);

} // namespace edgefem
