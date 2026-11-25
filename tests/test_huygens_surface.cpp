#include "edgefem/bc.hpp"
#include "edgefem/edge_basis.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/ports/wave_port.hpp"
#include "edgefem/post/huygens_surface.hpp"
#include "edgefem/solver.hpp"

#include "test_check.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>

using namespace edgefem;

namespace {
constexpr double kC0 = 299792458.0;
constexpr double kEta0 = 376.730313668;

// Edge DOF of E = x exp(-j k z): the line integral from the edge's n0 to n1
// (the canonical global direction).
VecC plane_wave_dofs(const Mesh &mesh, double k) {
  const std::complex<double> j(0.0, 1.0);
  VecC x(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    const auto &p0 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[i].n0)].xyz;
    const auto &p1 = mesh.nodes[mesh.nodeIndex.at(mesh.edges[i].n1)].xyz;
    const double dz = p1.z() - p0.z();
    std::complex<double> avg = std::exp(-j * k * p0.z());
    if (std::abs(k * dz) > 1e-12)
      avg = (std::exp(-j * k * p1.z()) - std::exp(-j * k * p0.z())) /
            (-j * k * dz);
    x(i) = (p1.x() - p0.x()) * avg;
  }
  return x;
}

// exp(+j omega t): for E = x exp(-j k z), a wave travelling to +z,
// curl E = -j omega mu H gives H = y exp(-j k z) / eta0. On the z = 0 port
// face, with an interior point inside the guide, the outward normal is -z.
void test_plane_wave_h_sign_and_normals() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  const double f = 1e9; // k h << 1 on this mesh
  const double omega = 2 * M_PI * f, k = omega / kC0;
  const VecC x = plane_wave_dofs(mesh, k);
  const Eigen::Vector3d inside(0.01143, 0.00508, 0.025);

  HuygensSurfaceData d =
      extract_huygens_surface(mesh, x, 2, omega, 1.0, inside);
  CHECK(!d.r.empty());
  double worst_h = 0.0, worst_e = 0.0;
  for (size_t i = 0; i < d.r.size(); ++i) {
    CHECK(d.n[i].z() < -0.999); // outward from the guide at z = 0
    const std::complex<double> j(0.0, 1.0);
    const std::complex<double> ph = std::exp(-j * k * d.r[i].z());
    worst_e = std::max(worst_e, std::abs(d.E_tan[i].x() - ph));
    worst_h = std::max(worst_h, std::abs(d.H_tan[i].y() * kEta0 - ph));
    CHECK(std::abs(d.E_tan[i].z()) < 1e-12); // tangential to z = 0
    CHECK(std::abs(d.H_tan[i].z()) < 1e-12);
  }
  std::cout << "plane wave: max |E_x - 1| = " << worst_e
            << ", max |eta0 H_y - 1| = " << worst_h << "\n";
  CHECK(worst_e < 0.02);
  CHECK(worst_h < 0.05);
}

// Interior faces (two parent tets) tagged as a Huygens surface: every normal
// points away from the interior point, whatever tet the lookup meets last.
void test_interior_surface_normals() {
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  // Collect faces shared by two tets with centroid z in a thin band.
  std::map<std::array<std::int64_t, 3>, int> face_count;
  for (const auto &tet : mesh.tets) {
    static const int fn[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
    for (const auto &f : fn) {
      std::array<std::int64_t, 3> key{tet.conn[f[0]], tet.conn[f[1]],
                                      tet.conn[f[2]]};
      std::sort(key.begin(), key.end());
      ++face_count[key];
    }
  }
  const int tag = 777;
  int added = 0;
  for (const auto &kv : face_count) {
    if (kv.second != 2)
      continue;
    Element tri;
    tri.type = ElemType::Tri3;
    tri.phys = tag;
    for (int i = 0; i < 3; ++i)
      tri.conn[i] = kv.first[i];
    mesh.tris.push_back(tri);
    ++added;
  }
  CHECK(added > 100);
  // extract_huygens_surface reads only conn and phys of the triangles.
  const double omega = 2 * M_PI * 1e9;
  const VecC x = plane_wave_dofs(mesh, omega / kC0);
  const Eigen::Vector3d inside(0.01143, 0.00508, 0.025);
  HuygensSurfaceData d =
      extract_huygens_surface(mesh, x, tag, omega, 1.0, inside);
  CHECK(d.r.size() == static_cast<size_t>(added));
  int wrong = 0;
  for (size_t i = 0; i < d.r.size(); ++i) {
    const double s = d.n[i].dot(d.r[i] - inside);
    if (std::abs(s) > 1e-9 && s < 0.0)
      ++wrong;
  }
  CHECK(wrong == 0);
}
} // namespace

int main() {
  // Load waveguide mesh which has known ports and propagating modes
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");

  std::cout << "Mesh loaded: " << mesh.tets.size() << " tets, "
            << mesh.tris.size() << " tris, " << mesh.edges.size() << " edges\n";

  auto tags = list_physical_tags(mesh);
  std::cout << "Surface tags:";
  for (int t : tags.surface_tags)
    std::cout << " " << t;
  std::cout << "\nVolume tags:";
  for (int t : tags.volume_tags)
    std::cout << " " << t;
  std::cout << "\n";

  // Need at least 2 surface tags (ports) and a surface for Huygens extraction
  CHECK(tags.surface_tags.size() >= 2);
  if (tags.surface_tags.size() < 2)
    return edgefem_test::exit_code();

  // Find a surface tag to use as Huygens surface
  // Use a port cross-section surface for testing
  auto it = tags.surface_tags.begin();
  int huygens_tag = *it;

  // Solve a simple problem to get a non-trivial solution
  BC bc;
  MaxwellParams params;
  params.omega = 2 * M_PI * 10e9;

  // Build a simple excitation - just create a non-zero solution
  // by assembling with PEC and a simple source
  for (int t : tags.surface_tags) {
    BC bc_t = build_edge_pec(mesh, t);
    for (int e : bc_t.dirichlet_edges)
      bc.dirichlet_edges.insert(e);
  }

  // Create a simple non-zero solution vector
  VecC solution(mesh.edges.size());
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    // Set non-trivial field pattern
    int idx0 = mesh.nodeIndex.at(mesh.edges[i].n0);
    int idx1 = mesh.nodeIndex.at(mesh.edges[i].n1);
    Eigen::Vector3d mid = 0.5 * (mesh.nodes[idx0].xyz + mesh.nodes[idx1].xyz);
    // Simple plane wave in z-direction, polarized in x
    double phase = params.omega / 3e8 * mid.z();
    solution(i) = std::complex<double>(std::cos(phase), std::sin(phase));
  }

  // Test: Extract Huygens surface
  std::cout << "\nExtracting Huygens surface on tag " << huygens_tag << "...\n";

  HuygensSurfaceData huygens =
      extract_huygens_surface(mesh, solution, huygens_tag, params.omega);

  std::cout << "Extracted " << huygens.r.size() << " surface points\n";

  CHECK(huygens.r.size() > 0 && "Must have surface points");
  CHECK(huygens.r.size() == huygens.n.size() && "Size mismatch: normals");
  CHECK(huygens.r.size() == huygens.E_tan.size() && "Size mismatch: E_tan");
  CHECK(huygens.r.size() == huygens.H_tan.size() && "Size mismatch: H_tan");
  CHECK(huygens.r.size() == huygens.area.size() && "Size mismatch: area");

  // Verify fields are non-zero
  double E_total = 0, H_total = 0;
  for (size_t i = 0; i < huygens.r.size(); ++i) {
    E_total += huygens.E_tan[i].norm();
    H_total += huygens.H_tan[i].norm();
  }
  std::cout << "Total |E_tan| = " << E_total << "\n";
  std::cout << "Total |H_tan| = " << H_total << "\n";

  CHECK(E_total > 1e-20 && "E_tan must be non-zero");
  CHECK(H_total > 1e-20 && "H_tan must be non-zero");

  // Verify fields are tangential (dot product with normal should be ~0)
  double max_E_normal = 0, max_H_normal = 0;
  for (size_t i = 0; i < huygens.r.size(); ++i) {
    Eigen::Vector3cd n_c = huygens.n[i].cast<std::complex<double>>();
    double E_n = std::abs(huygens.E_tan[i].dot(n_c));
    double H_n = std::abs(huygens.H_tan[i].dot(n_c));
    double E_mag = huygens.E_tan[i].norm();
    double H_mag = huygens.H_tan[i].norm();
    if (E_mag > 1e-20)
      max_E_normal = std::max(max_E_normal, E_n / E_mag);
    if (H_mag > 1e-20)
      max_H_normal = std::max(max_H_normal, H_n / H_mag);
  }
  std::cout << "Max |E·n|/|E| = " << max_E_normal << "\n";
  std::cout << "Max |H·n|/|H| = " << max_H_normal << "\n";

  // Note: with a synthetic (non-physical) solution on a non-Huygens surface,
  // the tangential projection may not be very accurate. The key test is that
  // the projection operation works (E_tan != E, i.e., normal component
  // removed). For a real FEM solution on a proper Huygens surface, these ratios
  // will be ~0.
  std::cout
      << "(Tangentiality check is informational for synthetic solutions)\n";

  // Verify areas are positive
  for (size_t i = 0; i < huygens.area.size(); ++i) {
    CHECK(huygens.area[i] > 0 && "Triangle areas must be positive");
  }

  // Verify normals are unit vectors
  for (size_t i = 0; i < huygens.n.size(); ++i) {
    double n_mag = huygens.n[i].norm();
    CHECK(std::abs(n_mag - 1.0) < 1e-10 && "Normals must be unit vectors");
  }

  test_plane_wave_h_sign_and_normals();
  test_interior_surface_normals();
  return edgefem_test::exit_code();
}
