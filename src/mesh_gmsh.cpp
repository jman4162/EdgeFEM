#include "vectorem/mesh.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace vectorem {

Mesh load_gmsh_v2(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open mesh file: " + path);
  }

  Mesh mesh;
  std::string line;
  while (std::getline(in, line)) {
    if (line == "$Nodes") {
      int n;
      in >> n;
      for (int i = 0; i < n; ++i) {
        std::int64_t id;
        double x, y, z;
        in >> id >> x >> y >> z;
        mesh.nodeIndex[id] = static_cast<int>(mesh.nodes.size());
        mesh.nodes.push_back(Node{id, Eigen::Vector3d{x, y, z}});
      }
      std::getline(in, line); // consume endline
      std::getline(in, line); // "$EndNodes"
    } else if (line == "$Elements") {
      int m;
      in >> m;
      for (int i = 0; i < m; ++i) {
        std::int64_t id;
        int type, num_tags;
        in >> id >> type >> num_tags;
        Element e;
        e.id = id;
        e.type = static_cast<ElemType>(type);
        for (int t = 0; t < num_tags; ++t) {
          int tag;
          in >> tag;
          if (t == 0)
            e.phys = tag;
        }
        if (type == 2) { // Tri3
          for (int k = 0; k < 3; ++k)
            in >> e.conn[k];
          mesh.tris.push_back(e);
        } else if (type == 4) { // Tet4
          for (int k = 0; k < 4; ++k)
            in >> e.conn[k];
          mesh.tets.push_back(e);
        } else {
          std::string dummy;
          std::getline(in, dummy);
        }
      }
      std::getline(in, line); // consume endline
      std::getline(in, line); // "$EndElements"
    }
  }
  return mesh;
}

} // namespace vectorem
