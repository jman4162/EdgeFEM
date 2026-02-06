#include <cassert>
#include <cmath>
#include <iostream>

#include "edgefem/mesh.hpp"
#include "edgefem/mesh_quality.hpp"

using namespace edgefem;

void test_tet_volume_positive() {
  std::cout << "test_tet_volume_positive..." << std::endl;

  // Regular tetrahedron with unit edge
  Eigen::Vector3d p0(0, 0, 0);
  Eigen::Vector3d p1(1, 0, 0);
  Eigen::Vector3d p2(0.5, std::sqrt(3.0) / 2.0, 0);
  Eigen::Vector3d p3(0.5, std::sqrt(3.0) / 6.0, std::sqrt(2.0 / 3.0));

  double vol = tet_volume(p0, p1, p2, p3);

  // Volume of regular tet with edge a is a^3 / (6*sqrt(2))
  double expected = 1.0 / (6.0 * std::sqrt(2.0));
  assert(std::abs(vol - expected) < 1e-10);

  std::cout << "  Volume: " << vol << " (expected: " << expected << ")" << std::endl;
}

void test_tet_volume_inverted() {
  std::cout << "test_tet_volume_inverted..." << std::endl;

  // Same tet but with inverted orientation
  Eigen::Vector3d p0(0, 0, 0);
  Eigen::Vector3d p1(1, 0, 0);
  Eigen::Vector3d p2(0.5, std::sqrt(3.0) / 2.0, 0);
  Eigen::Vector3d p3(0.5, std::sqrt(3.0) / 6.0, -std::sqrt(2.0 / 3.0));  // Negated z

  double vol = tet_volume(p0, p1, p2, p3);

  assert(vol < 0);
  std::cout << "  Inverted tet volume: " << vol << " (correctly negative)" << std::endl;
}

void test_tet_volume_degenerate() {
  std::cout << "test_tet_volume_degenerate..." << std::endl;

  // All points on a plane (degenerate)
  Eigen::Vector3d p0(0, 0, 0);
  Eigen::Vector3d p1(1, 0, 0);
  Eigen::Vector3d p2(0, 1, 0);
  Eigen::Vector3d p3(1, 1, 0);  // Coplanar

  double vol = tet_volume(p0, p1, p2, p3);

  assert(std::abs(vol) < 1e-15);
  std::cout << "  Degenerate tet volume: " << vol << " (correctly ~0)" << std::endl;
}

void test_validate_mesh_good() {
  std::cout << "test_validate_mesh_good..." << std::endl;

  // Load a known-good mesh
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");

  auto report = validate_mesh(mesh);

  std::cout << "  Tets: " << mesh.tets.size() << std::endl;
  std::cout << "  Degenerate: " << report.num_degenerate_tets << std::endl;
  std::cout << "  Inverted: " << report.num_inverted_tets << std::endl;
  std::cout << "  Min volume: " << report.min_volume << std::endl;
  std::cout << "  Max aspect ratio: " << report.max_aspect_ratio << std::endl;

  assert(report.is_valid());
  assert(report.min_volume > 0);
  assert(report.max_aspect_ratio > 0);
  assert(report.max_aspect_ratio < 100);  // Reasonable aspect ratio
}

void test_aspect_ratio() {
  std::cout << "test_aspect_ratio..." << std::endl;

  // Create a mesh with a single tet for testing
  Mesh mesh;

  // Add nodes
  mesh.nodes.push_back({0, {0, 0, 0}});
  mesh.nodes.push_back({1, {1, 0, 0}});
  mesh.nodes.push_back({2, {0, 1, 0}});
  mesh.nodes.push_back({3, {0, 0, 1}});

  mesh.nodeIndex[0] = 0;
  mesh.nodeIndex[1] = 1;
  mesh.nodeIndex[2] = 2;
  mesh.nodeIndex[3] = 3;

  Element tet;
  tet.id = 0;
  tet.type = ElemType::Tet4;
  tet.conn = {0, 1, 2, 3};
  mesh.tets.push_back(tet);

  auto report = validate_mesh(mesh);

  // For a regular tet, aspect ratio should be 1.0
  // For this right-angle tet, longest edge is sqrt(2), shortest is 1
  double expected_ar = std::sqrt(2.0);
  assert(std::abs(report.max_aspect_ratio - expected_ar) < 1e-10);

  std::cout << "  Aspect ratio: " << report.max_aspect_ratio << std::endl;
}

int main() {
  test_tet_volume_positive();
  test_tet_volume_inverted();
  test_tet_volume_degenerate();
  test_validate_mesh_good();
  test_aspect_ratio();
  std::cout << "All mesh quality tests passed!" << std::endl;
  return 0;
}
