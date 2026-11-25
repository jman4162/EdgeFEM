// load_gmsh_v2 input handling: it must read Gmsh 2.x ASCII files with LF or
// CRLF line endings, and refuse what it cannot read (MSH 4, binary MSH,
// second-order elements, truncated files) with an error instead of
// returning a wrong or empty mesh.

#include "test_check.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "edgefem/mesh.hpp"

using namespace edgefem;

namespace {
std::string read_file(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::string write_temp(const std::string &name, const std::string &text) {
  const std::string path = "/tmp/edgefem_test_" + name + ".msh";
  std::ofstream out(path, std::ios::binary);
  out << text;
  return path;
}

// Throws std::runtime_error whose message contains `needle`.
bool throws_with(const std::string &path, const std::string &needle) {
  try {
    (void)load_gmsh_v2(path);
  } catch (const std::runtime_error &e) {
    const bool ok = std::string(e.what()).find(needle) != std::string::npos;
    if (!ok)
      std::cout << "  message: " << e.what() << "\n";
    return ok;
  }
  return false;
}

std::string replace_all(std::string s, const std::string &from,
                        const std::string &to) {
  for (size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos;
       pos += to.size())
    s.replace(pos, from.size(), to);
  return s;
}

const char *kTinyV2 = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
                      "$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n"
                      "$EndNodes\n"
                      "$Elements\n2\n1 2 2 7 1 1 2 3\n2 4 2 9 1 1 2 3 4\n"
                      "$EndElements\n";
} // namespace

int main() {
  const Mesh ref = load_gmsh_v2("examples/cube_cavity.msh");
  CHECK(!ref.tets.empty());

  // CRLF line endings load the same mesh.
  {
    const std::string lf = read_file("examples/cube_cavity.msh");
    const std::string path = write_temp("crlf", replace_all(lf, "\n", "\r\n"));
    const Mesh m = load_gmsh_v2(path);
    CHECK(m.nodes.size() == ref.nodes.size());
    CHECK(m.tets.size() == ref.tets.size());
    CHECK(m.tris.size() == ref.tris.size());
    CHECK(m.edges.size() == ref.edges.size());
  }

  // A minimal valid file: one tet, one tri, physical tags read.
  {
    const Mesh m = load_gmsh_v2(write_temp("tiny", kTinyV2));
    CHECK(m.nodes.size() == 4 && m.tets.size() == 1 && m.tris.size() == 1);
    CHECK(m.tets.size() == 1 && m.tets[0].phys == 9);
    CHECK(m.tris.size() == 1 && m.tris[0].phys == 7);
  }

  // MSH 4.1 (Gmsh's default output) is refused with the msh2 hint.
  {
    const std::string v4 = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
                           "$Entities\n0 0 0 1\n1 0 0 0 1 1 1 0 0\n"
                           "$EndEntities\n$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n"
                           "4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n";
    CHECK(throws_with(write_temp("v4", v4), "-format msh2"));
  }

  // Binary MSH 2 is refused.
  {
    const std::string bin = replace_all(kTinyV2, "2.2 0 8", "2.2 1 8");
    CHECK(throws_with(write_temp("binary", bin), "binary"));
  }

  // No $MeshFormat section at all.
  {
    const std::string none =
        std::string(kTinyV2).substr(std::string(kTinyV2).find("$Nodes"));
    CHECK(throws_with(write_temp("noformat", none), "$MeshFormat"));
  }

  // Second-order tets (type 11) are refused instead of dropped.
  {
    std::string t10 =
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n10\n"
        "1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n5 .5 0 0\n6 .5 .5 0\n7 0 .5 0\n"
        "8 0 0 .5\n9 0 .5 .5\n10 .5 0 .5\n$EndNodes\n"
        "$Elements\n1\n1 11 2 9 1 1 2 3 4 5 6 7 8 9 10\n$EndElements\n";
    CHECK(throws_with(write_temp("tet10", t10), "element type 11"));
  }

  // Truncated $Nodes block.
  {
    const std::string cut = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
                            "$Nodes\n4\n1 0 0 0\n2 1 0\n";
    CHECK(throws_with(write_temp("truncated", cut), "$Nodes"));
  }

  return edgefem_test::exit_code();
}
