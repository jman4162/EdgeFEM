#include "edgefem/bc.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace edgefem {

BC build_scalar_pec(const Mesh &mesh, int pec_tag) {
  BC bc;
  for (const auto &tri : mesh.tris) {
    if (tri.phys == pec_tag) {
      for (int k = 0; k < 3; ++k) {
        auto it = mesh.nodeIndex.find(tri.conn[k]);
        if (it != mesh.nodeIndex.end()) {
          bc.dirichlet_nodes.insert(it->second);
        }
      }
    }
  }

  // Warn if no nodes found (common user error: wrong tag)
  if (bc.dirichlet_nodes.empty()) {
    auto tags = list_physical_tags(mesh);
    std::ostringstream oss;
    oss << "WARNING: build_scalar_pec() found no nodes for PEC tag " << pec_tag << ".\n";
    oss << "  Available surface tags: ";
    if (tags.surface_tags.empty()) {
      oss << "(none)";
    } else {
      bool first = true;
      for (int t : tags.surface_tags) {
        if (!first) oss << ", ";
        oss << t;
        first = false;
      }
    }
    oss << "\n";
    std::cerr << oss.str();
  }

  return bc;
}

BC build_edge_pec(const Mesh &mesh, int pec_tag) {
  BC bc;
  for (const auto &tri : mesh.tris) {
    if (tri.phys == pec_tag) {
      for (int k = 0; k < 3; ++k) {
        bc.dirichlet_edges.insert(tri.edges[k]);
      }
    }
  }

  // Warn if no edges found (common user error: wrong tag)
  if (bc.dirichlet_edges.empty()) {
    auto tags = list_physical_tags(mesh);
    std::ostringstream oss;
    oss << "WARNING: build_edge_pec() found no edges for PEC tag " << pec_tag << ".\n";
    oss << "  Available surface tags: ";
    if (tags.surface_tags.empty()) {
      oss << "(none)";
    } else {
      bool first = true;
      for (int t : tags.surface_tags) {
        if (!first) oss << ", ";
        oss << t;
        first = false;
      }
    }
    oss << "\n";
    std::cerr << oss.str();
  }

  return bc;
}

PhysicalTagInfo list_physical_tags(const Mesh &mesh) {
  PhysicalTagInfo info;
  for (const auto &tet : mesh.tets) {
    info.volume_tags.insert(tet.phys);
  }
  for (const auto &tri : mesh.tris) {
    info.surface_tags.insert(tri.phys);
  }
  return info;
}

bool has_surface_tag(const Mesh &mesh, int tag) {
  for (const auto &tri : mesh.tris) {
    if (tri.phys == tag) return true;
  }
  return false;
}

bool has_volume_tag(const Mesh &mesh, int tag) {
  for (const auto &tet : mesh.tets) {
    if (tet.phys == tag) return true;
  }
  return false;
}

} // namespace edgefem
