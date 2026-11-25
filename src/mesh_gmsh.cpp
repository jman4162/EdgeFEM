#include "edgefem/mesh.hpp"
#include "edgefem/mesh_quality.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace edgefem {
namespace {
void build_edges(Mesh &mesh);
}

namespace {
// getline that drops a trailing '\r' (CRLF files).
bool get_line(std::istream &in, std::string &line) {
  if (!std::getline(in, line))
    return false;
  if (!line.empty() && line.back() == '\r')
    line.pop_back();
  return true;
}

[[noreturn]] void fail(const std::string &path, const std::string &what) {
  throw std::runtime_error("load_gmsh_v2: " + path + ": " + what);
}

// Skips lines up to and including `end_marker`.
void skip_to(std::istream &in, const std::string &end_marker,
             const std::string &path) {
  std::string line;
  while (get_line(in, line))
    if (line == end_marker)
      return;
  fail(path, "missing " + end_marker);
}
} // namespace

Mesh load_gmsh_v2(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open mesh file: " + path);
  }

  Mesh mesh;
  std::string line;
  bool have_format = false;
  while (get_line(in, line)) {
    if (line == "$MeshFormat") {
      if (!get_line(in, line))
        fail(path, "truncated $MeshFormat");
      std::istringstream fs(line);
      double version = 0.0;
      int file_type = -1, data_size = 0;
      if (!(fs >> version >> file_type >> data_size))
        fail(path, "unreadable $MeshFormat line '" + line + "'");
      if (version < 2.0 || version >= 3.0)
        fail(path, "MSH version " + line.substr(0, line.find(' ')) +
                       " is not supported; only MSH 2.x ASCII is. Export "
                       "with `gmsh -format msh2` (or Mesh.MshFileVersion = "
                       "2.2)");
      if (file_type != 0)
        fail(path, "binary MSH is not supported; export ASCII with "
                   "`gmsh -format msh2` and without -bin");
      skip_to(in, "$EndMeshFormat", path);
      have_format = true;
    } else if (line == "$Nodes") {
      if (!have_format)
        fail(path, "no $MeshFormat section before $Nodes (not an MSH 2 "
                   "file?)");
      if (!get_line(in, line))
        fail(path, "truncated $Nodes");
      long n = -1;
      std::istringstream(line) >> n;
      if (n < 0)
        fail(path, "bad node count in $Nodes: '" + line + "'");
      mesh.nodes.reserve(static_cast<size_t>(n));
      for (long i = 0; i < n; ++i) {
        std::int64_t id;
        double x, y, z;
        if (!get_line(in, line) ||
            !(std::istringstream(line) >> id >> x >> y >> z))
          fail(path, "truncated or malformed $Nodes at node " +
                         std::to_string(i + 1) + " of " + std::to_string(n));
        mesh.nodeIndex[id] = static_cast<int>(mesh.nodes.size());
        mesh.nodes.push_back(Node{id, Eigen::Vector3d{x, y, z}});
      }
      if (!get_line(in, line) || line != "$EndNodes")
        fail(path, "expected $EndNodes after " + std::to_string(n) + " nodes");
    } else if (line == "$Elements") {
      if (!have_format)
        fail(path, "no $MeshFormat section before $Elements");
      if (!get_line(in, line))
        fail(path, "truncated $Elements");
      long m = -1;
      std::istringstream(line) >> m;
      if (m < 0)
        fail(path, "bad element count in $Elements: '" + line + "'");
      for (long i = 0; i < m; ++i) {
        if (!get_line(in, line))
          fail(path, "truncated $Elements at element " + std::to_string(i + 1) +
                         " of " + std::to_string(m));
        std::istringstream es(line);
        std::int64_t id;
        int type, num_tags;
        if (!(es >> id >> type >> num_tags) || num_tags < 0)
          fail(path, "malformed element line '" + line + "'");
        Element e;
        e.id = id;
        e.type = static_cast<ElemType>(type);
        for (int t = 0; t < num_tags; ++t) {
          int tag;
          if (!(es >> tag))
            fail(path, "malformed element line '" + line + "'");
          if (t == 0)
            e.phys = tag;
        }
        bool ok = true;
        if (type == 2) { // Tri3
          for (int k = 0; k < 3; ++k)
            ok = ok && static_cast<bool>(es >> e.conn[k]);
          mesh.tris.push_back(e);
        } else if (type == 1) { // Line2
          BoundaryLine b;
          b.phys = e.phys;
          ok = static_cast<bool>(es >> b.n0 >> b.n1);
          mesh.boundary_lines.push_back(b);
        } else if (type == 4) { // Tet4
          for (int k = 0; k < 4; ++k)
            ok = ok && static_cast<bool>(es >> e.conn[k]);
          mesh.tets.push_back(e);
        } else if (type == 15) { // point: not used
        } else {
          fail(path, "element type " + std::to_string(type) +
                         " is not supported (Line2 = 1, Tri3 = 2, Tet4 = 4 "
                         "and points only; mesh with -order 1)");
        }
        if (!ok)
          fail(path, "malformed element line '" + line + "'");
      }
      if (!get_line(in, line) || line != "$EndElements")
        fail(path,
             "expected $EndElements after " + std::to_string(m) + " elements");
    }
  }
  if (!have_format)
    fail(path, "no $MeshFormat section (not an MSH 2 file)");
  for (const auto &t : mesh.tets)
    for (int k = 0; k < 4; ++k)
      if (!mesh.nodeIndex.count(t.conn[k]))
        fail(path, "tet " + std::to_string(t.id) + " uses undefined node " +
                       std::to_string(t.conn[k]));
  for (const auto &t : mesh.tris)
    for (int k = 0; k < 3; ++k)
      if (!mesh.nodeIndex.count(t.conn[k]))
        fail(path, "triangle " + std::to_string(t.id) +
                       " uses undefined node " + std::to_string(t.conn[k]));
  build_edges(mesh);

  // Validate mesh quality and warn about issues
  if (!mesh.tets.empty()) {
    auto report = validate_mesh(mesh);
    if (report.num_degenerate_tets > 0) {
      std::cerr << "WARNING: Mesh '" << path << "' contains "
                << report.num_degenerate_tets
                << " degenerate (zero-volume) tetrahedra. "
                   "Results may be inaccurate."
                << std::endl;
    }
    if (report.num_inverted_tets > 0) {
      std::cerr << "WARNING: Mesh '" << path << "' contains "
                << report.num_inverted_tets
                << " inverted tetrahedra. Results may be inaccurate."
                << std::endl;
    }
    if (report.max_aspect_ratio > 100.0) {
      std::cerr << "WARNING: Mesh '" << path << "' has extreme aspect ratio ("
                << report.max_aspect_ratio
                << "). Consider refining the mesh for better accuracy."
                << std::endl;
    }
  }

  return mesh;
}

namespace {
void build_edges(Mesh &mesh) {
  const std::array<std::array<int, 2>, 6> tet_pairs = {
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2},
      std::array<int, 2>{0, 3}, std::array<int, 2>{1, 2},
      std::array<int, 2>{1, 3}, std::array<int, 2>{2, 3}};
  const std::array<std::array<int, 2>, 3> tri_pairs = {
      std::array<int, 2>{0, 1}, std::array<int, 2>{1, 2},
      std::array<int, 2>{2, 0}};
  for (auto &tet : mesh.tets) {
    for (int e = 0; e < 6; ++e) {
      auto a = tet.conn[tet_pairs[e][0]];
      auto b = tet.conn[tet_pairs[e][1]];
      int sign = (a < b) ? 1 : -1;
      std::uint64_t key = make_edge_key(a, b);
      auto it = mesh.edgeIndex.find(key);
      if (it == mesh.edgeIndex.end()) {
        int idx = static_cast<int>(mesh.edges.size());
        mesh.edgeIndex[key] = idx;
        mesh.edges.push_back(Edge{std::min(a, b), std::max(a, b)});
        it = mesh.edgeIndex.find(key);
      }
      tet.edges[e] = it->second;
      tet.edge_orient[e] = sign;
    }
  }
  for (auto &tri : mesh.tris) {
    for (int e = 0; e < 3; ++e) {
      auto a = tri.conn[tri_pairs[e][0]];
      auto b = tri.conn[tri_pairs[e][1]];
      int sign = (a < b) ? 1 : -1;
      std::uint64_t key = make_edge_key(a, b);
      auto it = mesh.edgeIndex.find(key);
      if (it == mesh.edgeIndex.end()) {
        int idx = static_cast<int>(mesh.edges.size());
        mesh.edgeIndex[key] = idx;
        mesh.edges.push_back(Edge{std::min(a, b), std::max(a, b)});
        it = mesh.edgeIndex.find(key);
      }
      tri.edges[e] = it->second;
      tri.edge_orient[e] = sign;
    }
  }
}
} // namespace
} // namespace edgefem
