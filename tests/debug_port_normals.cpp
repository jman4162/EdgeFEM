// Debug port normal directions and weight signs
#include "edgefem/maxwell.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace edgefem;

int main() {
  std::cout << std::setprecision(8);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");

  // Check triangle normals for both ports
  std::cout << "=== Port 1 (tag 2, z=0) Triangle Normals ===" << std::endl;
  double port1_normal_sum = 0;
  int port1_count = 0;
  for (const auto &tri : mesh.tris) {
    if (tri.phys != 2)
      continue;

    const auto &n0 = mesh.nodes[mesh.nodeIndex.at(tri.conn[0])].xyz;
    const auto &n1 = mesh.nodes[mesh.nodeIndex.at(tri.conn[1])].xyz;
    const auto &n2 = mesh.nodes[mesh.nodeIndex.at(tri.conn[2])].xyz;

    Eigen::Vector3d v1 = n1 - n0;
    Eigen::Vector3d v2 = n2 - n0;
    Eigen::Vector3d normal = v1.cross(v2);

    port1_normal_sum += normal.z();
    port1_count++;

    if (port1_count <= 5) {
      std::cout << "  Tri " << port1_count << ": z=" << n0.z() << ", normal=("
                << normal.x() << "," << normal.y() << "," << normal.z() << ")"
                << std::endl;
    }
  }
  std::cout << "  Total triangles: " << port1_count << std::endl;
  std::cout << "  Sum of normal_z: " << port1_normal_sum << std::endl;
  std::cout << "  Port 1 normal sign: " << (port1_normal_sum >= 0 ? "+1" : "-1")
            << std::endl;

  std::cout << "\n=== Port 2 (tag 3, z=0.05) Triangle Normals ===" << std::endl;
  double port2_normal_sum = 0;
  int port2_count = 0;
  for (const auto &tri : mesh.tris) {
    if (tri.phys != 3)
      continue;

    const auto &n0 = mesh.nodes[mesh.nodeIndex.at(tri.conn[0])].xyz;
    const auto &n1 = mesh.nodes[mesh.nodeIndex.at(tri.conn[1])].xyz;
    const auto &n2 = mesh.nodes[mesh.nodeIndex.at(tri.conn[2])].xyz;

    Eigen::Vector3d v1 = n1 - n0;
    Eigen::Vector3d v2 = n2 - n0;
    Eigen::Vector3d normal = v1.cross(v2);

    port2_normal_sum += normal.z();
    port2_count++;

    if (port2_count <= 5) {
      std::cout << "  Tri " << port2_count << ": z=" << n0.z() << ", normal=("
                << normal.x() << "," << normal.y() << "," << normal.z() << ")"
                << std::endl;
    }
  }
  std::cout << "  Total triangles: " << port2_count << std::endl;
  std::cout << "  Sum of normal_z: " << port2_normal_sum << std::endl;
  std::cout << "  Port 2 normal sign: " << (port2_normal_sum >= 0 ? "+1" : "-1")
            << std::endl;

  // Build wave ports and check weight signs
  double freq = 10e9;
  RectWaveguidePort dims{0.02286, 0.01016};

  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  // Check weight statistics
  std::cout << "\n=== Port 1 Weight Statistics ===" << std::endl;
  double w1_real_sum = 0, w1_imag_sum = 0;
  double w1_real_pos = 0, w1_real_neg = 0;
  double w1_imag_pos = 0, w1_imag_neg = 0;
  for (int i = 0; i < wp1.weights.size(); ++i) {
    double re = std::real(wp1.weights(i));
    double im = std::imag(wp1.weights(i));
    w1_real_sum += re;
    w1_imag_sum += im;
    if (re > 0)
      w1_real_pos += re;
    else
      w1_real_neg += re;
    if (im > 0)
      w1_imag_pos += im;
    else
      w1_imag_neg += im;
  }
  std::cout << "  Sum of real parts: " << w1_real_sum << std::endl;
  std::cout << "  Sum of imag parts: " << w1_imag_sum << std::endl;
  std::cout << "  Positive imag sum: " << w1_imag_pos << std::endl;
  std::cout << "  Negative imag sum: " << w1_imag_neg << std::endl;

  std::cout << "\n=== Port 2 Weight Statistics ===" << std::endl;
  double w2_real_sum = 0, w2_imag_sum = 0;
  double w2_imag_pos = 0, w2_imag_neg = 0;
  for (int i = 0; i < wp2.weights.size(); ++i) {
    double re = std::real(wp2.weights(i));
    double im = std::imag(wp2.weights(i));
    w2_real_sum += re;
    w2_imag_sum += im;
    if (im > 0)
      w2_imag_pos += im;
    else
      w2_imag_neg += im;
  }
  std::cout << "  Sum of real parts: " << w2_real_sum << std::endl;
  std::cout << "  Sum of imag parts: " << w2_imag_sum << std::endl;
  std::cout << "  Positive imag sum: " << w2_imag_pos << std::endl;
  std::cout << "  Negative imag sum: " << w2_imag_neg << std::endl;

  // Compare weights - they should have similar magnitude patterns
  std::cout << "\n=== Weight Comparison ===" << std::endl;
  std::cout << "||w1||² = " << wp1.weights.squaredNorm() << std::endl;
  std::cout << "||w2||² = " << wp2.weights.squaredNorm() << std::endl;

  // Check if port normals should be opposite
  std::cout << "\n=== Physical Interpretation ===" << std::endl;
  std::cout << "For waveguide propagation in +z direction:" << std::endl;
  std::cout << "  Port 1 (z=0): outward normal should be -z (into domain)"
            << std::endl;
  std::cout << "  Port 2 (z=L): outward normal should be +z (out of domain)"
            << std::endl;
  std::cout
      << "  But for S-parameter extraction, both ports should use consistent"
      << std::endl;
  std::cout << "  'forward wave' direction (from port 1 to port 2)."
            << std::endl;

  return 0;
}
