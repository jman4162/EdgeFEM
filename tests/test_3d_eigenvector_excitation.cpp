// Test using the 3D FEM eigenvector directly as excitation source
// Hypothesis: if we excite with a source that matches the 3D eigenvector,
// the wave should propagate correctly

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double kc = M_PI / 0.02286;
  double kc_sq = kc * kc;
  double beta = std::sqrt(k0*k0 - kc_sq);

  // Compute 3D FEM eigenvector
  MaxwellParams p_low;
  p_low.omega = 1e-10;
  auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);

  double f_test = 1e9;
  double k0_test = 2*M_PI*f_test/c0;
  MaxwellParams p_test;
  p_test.omega = 2*M_PI*f_test;
  auto asmbl_test = assemble_maxwell(mesh, p_test, bc, {}, -1);

  int n = asmbl_low.A.rows();
  std::vector<int> free_to_orig;
  for (int i = 0; i < n; ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      free_to_orig.push_back(i);
    }
  }
  int n_free = free_to_orig.size();

  Eigen::MatrixXd K_free(n_free, n_free), M_free(n_free, n_free);
  for (int i = 0; i < n_free; ++i) {
    for (int j = 0; j < n_free; ++j) {
      int oi = free_to_orig[i], oj = free_to_orig[j];
      K_free(i, j) = std::real(asmbl_low.A.coeff(oi, oj));
      double A_test_ij = std::real(asmbl_test.A.coeff(oi, oj));
      M_free(i, j) = (K_free(i, j) - A_test_ij) / (k0_test * k0_test);
    }
  }

  std::cout << "Computing 3D eigenvalue problem..." << std::endl;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // Find TE10
  int te10_idx = -1;
  double min_diff = 1e100;
  for (int i = 0; i < n_free; ++i) {
    if (eigenvalues(i) > 100) {
      double diff = std::abs(eigenvalues(i) - kc_sq);
      if (diff < min_diff) { min_diff = diff; te10_idx = i; }
    }
  }

  std::cout << "TE10 eigenvalue: " << eigenvalues(te10_idx) << " (expected: " << kc_sq << ")" << std::endl;

  // Create full eigenvector
  Eigen::VectorXd v_te10_real = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n_free; ++i) {
    v_te10_real(free_to_orig[i]) = eigenvectors.col(te10_idx)(i);
  }

  // Convert to complex
  Eigen::VectorXcd v_te10 = v_te10_real.cast<std::complex<double>>();

  // Normalize so ||v||² = 1
  v_te10 /= v_te10.norm();

  std::cout << "||v_te10|| = " << v_te10.norm() << std::endl;

  // Now set up the driven problem with eigenvector-based source
  // Use the operating frequency FEM matrix
  MaxwellParams p;
  p.omega = omega;

  // Assemble without ports (no loading)
  auto asmbl_driven = assemble_maxwell(mesh, p, bc, {}, -1);

  // Add ABC damping at all non-PEC edges on boundary surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  RectWaveguidePort dims{0.02286, 0.01016};
  PortMode mode1 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);

  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port2_surf, dims, mode2);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  std::complex<double> j(0, 1);
  std::complex<double> abc_damping = j * beta;

  // Add ABC to port edges
  for (int e : wp1.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      asmbl_driven.A.coeffRef(e, e) += abc_damping;
    }
  }
  for (int e : wp2.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      asmbl_driven.A.coeffRef(e, e) += abc_damping;
    }
  }

  // Create source from eigenvector restricted to port 1 edges
  // Scale so total source power is 1 W
  Eigen::VectorXcd b = Eigen::VectorXcd::Zero(n);

  double v_port_norm = 0;
  for (int e : wp1.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      v_port_norm += std::norm(v_te10(e));
    }
  }
  v_port_norm = std::sqrt(v_port_norm);

  double source_scale = 1.0 / v_port_norm;  // Normalize port source
  for (int e : wp1.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      b(e) = source_scale * v_te10(e);
    }
  }

  std::cout << "||b|| = " << b.norm() << std::endl;

  // Solve
  auto res = solve_linear(asmbl_driven.A, b, {});

  // Examine solution at ports
  std::cout << "\n=== Solution Analysis ===" << std::endl;

  // Solution at port 1
  std::complex<double> proj1 = 0;
  for (int e : wp1.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      proj1 += std::conj(v_te10(e)) * res.x(e);
    }
  }
  std::cout << "Projection at port 1: v^H * x = " << proj1 << std::endl;

  // Solution at port 2
  std::complex<double> proj2 = 0;
  for (int e : wp2.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      proj2 += std::conj(v_te10(e)) * res.x(e);
    }
  }
  std::cout << "Projection at port 2: v^H * x = " << proj2 << std::endl;

  // For a traveling wave, proj2/proj1 should be exp(-jβL)
  std::complex<double> ratio = proj2 / proj1;
  double L = 0.05;
  std::complex<double> expected_ratio = std::exp(std::complex<double>(0, -beta * L));

  std::cout << "\nRatio proj2/proj1 = " << ratio << " (|ratio| = " << std::abs(ratio) << ")" << std::endl;
  std::cout << "Expected ratio = " << expected_ratio << " (|·| = 1)" << std::endl;

  // Also check using analytical weights
  std::cout << "\n=== Using Analytical Weights ===" << std::endl;

  std::complex<double> V1 = 0, V2 = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V1 += std::conj(wp1.weights(i)) * res.x(e);
    }
  }
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      V2 += std::conj(wp2.weights(i)) * res.x(e);
    }
  }

  std::cout << "V1 (analytical weights) = " << V1 << std::endl;
  std::cout << "V2 (analytical weights) = " << V2 << std::endl;
  std::cout << "V2/V1 = " << V2/V1 << std::endl;

  // Examine solution magnitude along z
  std::cout << "\n=== Solution Magnitude vs Z ===" << std::endl;

  std::map<int, std::pair<double, int>> z_sum;  // z_mm -> (sum|x|, count)
  for (int e = 0; e < n; ++e) {
    if (bc.dirichlet_edges.count(e)) continue;

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    int z_mm = (int)((p0.z() + p1.z()) / 2.0 * 1000);

    z_sum[z_mm].first += std::abs(res.x(e));
    z_sum[z_mm].second++;
  }

  std::cout << "z(mm)\tmean|x|" << std::endl;
  for (const auto &kv : z_sum) {
    if (kv.first % 10 == 0 || kv.first == 0 || kv.first == 50) {
      std::cout << kv.first << "\t" << kv.second.first / kv.second.second << std::endl;
    }
  }

  return 0;
}
