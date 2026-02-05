// Test direct eigenmode excitation for S-parameter extraction
// Instead of using port weights, directly excite the 3D FEM eigenvector
// and compute S-parameters from the overlap integral

#include <iostream>
#include <iomanip>
#include <cmath>
#include <unordered_set>
#include "vectorem/maxwell.hpp"
#include "vectorem/solver.hpp"

using namespace vectorem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

// Compute the overlap integral of two edge fields on a port surface
// Returns ∫∫_S (E1 × H2*) · n dS which represents power flow
std::complex<double> compute_power_overlap(
    const Mesh &mesh,
    const PortSurfaceMesh &port_surf,
    const Eigen::VectorXcd &E1,
    const Eigen::VectorXcd &E2,
    const BC &bc,
    double omega,
    double mu) {

  std::complex<double> overlap = 0.0;

  // For edge elements, the overlap integral simplifies to
  // a sum over port edges weighted by edge lengths
  for (size_t tri_idx = 0; tri_idx < port_surf.mesh.tris.size(); ++tri_idx) {
    const auto &tri3d = mesh.tris.at(port_surf.volume_tri_indices.at(tri_idx));

    for (int e = 0; e < 3; ++e) {
      int edge_idx = tri3d.edges[e];
      if (bc.dirichlet_edges.count(edge_idx)) continue;

      // Get edge length
      const auto &edge = mesh.edges.at(edge_idx);
      const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
      double edge_len = (p1 - p0).norm();

      // Simplified overlap: E1* · E2 weighted by edge length
      overlap += std::conj(E1(edge_idx)) * E2(edge_idx) * edge_len;
    }
  }

  return overlap;
}

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);
  int m = mesh.edges.size();

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double L = 0.05;

  RectWaveguidePort dims{0.02286, 0.01016};
  double kc = M_PI / dims.a;
  double beta = std::sqrt(k0*k0 - kc*kc);
  double Z0_te10 = omega * mu0 / beta;  // TE mode impedance

  std::cout << "=== Direct Eigenmode Excitation Test ===" << std::endl;
  std::cout << "Frequency: " << freq/1e9 << " GHz" << std::endl;
  std::cout << "Waveguide length: " << L*1000 << " mm" << std::endl;
  std::cout << "TE10 beta = " << beta << " rad/m" << std::endl;
  std::cout << "TE10 Z0 = " << Z0_te10 << " ohms" << std::endl;
  std::cout << "Expected phase shift: " << -beta*L*180/M_PI << " deg" << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Collect port 1 and port 2 edges
  std::unordered_set<int> port1_edges, port2_edges;
  for (size_t tri_idx = 0; tri_idx < port1_surf.mesh.tris.size(); ++tri_idx) {
    const auto &tri3d = mesh.tris.at(port1_surf.volume_tri_indices.at(tri_idx));
    for (int e = 0; e < 3; ++e) {
      port1_edges.insert(tri3d.edges[e]);
    }
  }
  for (size_t tri_idx = 0; tri_idx < port2_surf.mesh.tris.size(); ++tri_idx) {
    const auto &tri3d = mesh.tris.at(port2_surf.volume_tri_indices.at(tri_idx));
    for (int e = 0; e < 3; ++e) {
      port2_edges.insert(tri3d.edges[e]);
    }
  }

  std::cout << "Port 1 edges: " << port1_edges.size() << std::endl;
  std::cout << "Port 2 edges: " << port2_edges.size() << std::endl;

  // Compute 3D FEM eigenvector for TE10 mode
  double kc_sq = kc * kc;
  Eigen::VectorXd v_te10 = compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);

  std::cout << "\nTE10 eigenvector norm: " << v_te10.norm() << std::endl;

  // Extract eigenvector at port surfaces
  Eigen::VectorXcd v1 = Eigen::VectorXcd::Zero(m);
  Eigen::VectorXcd v2 = Eigen::VectorXcd::Zero(m);

  for (int e : port1_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      v1(e) = std::complex<double>(0.0, v_te10(e));  // Multiply by j for TE mode
    }
  }
  for (int e : port2_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      v2(e) = std::complex<double>(0.0, v_te10(e));
    }
  }

  double v1_norm = v1.norm();
  double v2_norm = v2.norm();
  std::cout << "||v1|| = " << v1_norm << ", ||v2|| = " << v2_norm << std::endl;

  // Normalize eigenvector for unit amplitude
  if (v1_norm > 1e-10) v1 /= v1_norm;
  if (v2_norm > 1e-10) v2 /= v2_norm;

  // Assemble Maxwell system without ports
  MaxwellParams p;
  p.omega = omega;
  p.use_abc = false;
  p.use_port_abc = false;

  std::vector<WavePort> no_ports;
  auto asmbl = assemble_maxwell(mesh, p, bc, no_ports, -1);

  // Add ABC on port 1 and port 2 edges
  // The ABC for port edge e is: A(e,e) += j*beta * M_edge
  // where M_edge is the edge mass matrix contribution
  // For simplicity, use diagonal approximation: A(e,e) += j*beta
  const std::complex<double> abc_coeff(0.0, beta);

  for (int e : port1_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      asmbl.A.coeffRef(e, e) += abc_coeff;
    }
  }
  for (int e : port2_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      asmbl.A.coeffRef(e, e) += abc_coeff;
    }
  }

  // Create excitation: RHS = j*beta * v1 (incident mode at port 1)
  // The factor of 2 comes from wave equation with incident + reflected
  asmbl.b = 2.0 * abc_coeff * v1;

  // Solve
  auto res = solve_linear(asmbl.A, asmbl.b, {});
  std::cout << "\nSolver: " << res.method << ", iters=" << res.iters
            << ", converged=" << res.converged << std::endl;
  std::cout << "||x|| = " << res.x.norm() << std::endl;

  // Compute S-parameters from overlap integrals
  // S11 = <v1, E - E_inc> / <v1, E_inc> = <v1, E>/<v1, E_inc> - 1
  // S21 = <v2, E> / <v1, E_inc>

  // Port voltages via eigenvector overlap
  std::complex<double> V1 = v1.dot(res.x);
  std::complex<double> V2 = v2.dot(res.x);
  std::complex<double> V_inc = v1.dot(v1);  // = 1 since normalized

  std::cout << "\nV1 = " << V1 << std::endl;
  std::cout << "V2 = " << V2 << std::endl;
  std::cout << "V_inc = " << V_inc << std::endl;

  // S-parameters
  std::complex<double> S11 = (V1 - V_inc) / V_inc;
  std::complex<double> S21 = V2 / V_inc;

  std::cout << "\n=== S-Parameters ===" << std::endl;
  std::cout << "S11 = " << S11 << " (|S11| = " << std::abs(S11) << ")" << std::endl;
  std::cout << "S21 = " << S21 << " (|S21| = " << std::abs(S21)
            << ", phase = " << std::arg(S21)*180/M_PI << "°)" << std::endl;

  // Expected
  std::complex<double> S21_expected = std::exp(std::complex<double>(0, -beta * L));
  std::cout << "\nExpected:" << std::endl;
  std::cout << "S11 = 0" << std::endl;
  std::cout << "S21 = " << S21_expected << " (|S21| = 1, phase = "
            << std::arg(S21_expected)*180/M_PI << "°)" << std::endl;

  // Passivity check
  double passivity = std::norm(S11) + std::norm(S21);
  std::cout << "\nPassivity: |S11|² + |S21|² = " << passivity
            << " (should be ≤ 1)" << std::endl;

  // Check field at different z positions
  std::cout << "\n=== Field Along Waveguide ===" << std::endl;

  // Find edges at different z positions
  double z_min = 1e10, z_max = -1e10;
  for (const auto &node : mesh.nodes) {
    z_min = std::min(z_min, node.xyz.z());
    z_max = std::max(z_max, node.xyz.z());
  }

  for (double z_frac : {0.0, 0.25, 0.5, 0.75, 1.0}) {
    double z_target = z_min + z_frac * (z_max - z_min);

    // Find edges near this z
    std::complex<double> E_sum = 0.0;
    int count = 0;
    for (size_t i = 0; i < mesh.edges.size(); ++i) {
      if (bc.dirichlet_edges.count(i)) continue;

      const auto &edge = mesh.edges[i];
      const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
      double z_edge = (p0.z() + p1.z()) / 2.0;

      if (std::abs(z_edge - z_target) < 0.005) {
        E_sum += res.x(i);
        count++;
      }
    }

    if (count > 0) {
      E_sum /= count;
      std::cout << "z=" << z_target*1000 << "mm: avg |E| = " << std::abs(E_sum)
                << ", phase = " << std::arg(E_sum)*180/M_PI << "°" << std::endl;
    }
  }

  return 0;
}
