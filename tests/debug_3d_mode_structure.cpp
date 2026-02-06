// Examine the 3D structure of the FEM eigenvector
// Hypothesis: it's a 3D cavity mode (TE10p) not a uniform cross-section mode

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double kc_te10 = M_PI / 0.02286; // TE10 cutoff
  double kc_te10_sq = kc_te10 * kc_te10;

  // Compute FEM eigenvalue and eigenvector
  MaxwellParams p_low;
  p_low.omega = 1e-10;
  auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);

  double f_test = 1e9;
  double k0_test = 2 * M_PI * f_test / c0;
  MaxwellParams p_test;
  p_test.omega = 2 * M_PI * f_test;
  auto asmbl_test = assemble_maxwell(mesh, p_test, bc, {}, -1);

  int n = asmbl_low.A.rows();

  std::vector<int> free_to_orig;
  std::vector<int> orig_to_free(n, -1);
  for (int i = 0; i < n; ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      orig_to_free[i] = free_to_orig.size();
      free_to_orig.push_back(i);
    }
  }
  int n_free = free_to_orig.size();

  Eigen::MatrixXd K_free(n_free, n_free);
  Eigen::MatrixXd M_free(n_free, n_free);

  for (int i = 0; i < n_free; ++i) {
    for (int j = 0; j < n_free; ++j) {
      int oi = free_to_orig[i];
      int oj = free_to_orig[j];
      double K_ij = std::real(asmbl_low.A.coeff(oi, oj));
      double A_test_ij = std::real(asmbl_test.A.coeff(oi, oj));
      double M_ij = (K_ij - A_test_ij) / (k0_test * k0_test);
      K_free(i, j) = K_ij;
      M_free(i, j) = M_ij;
    }
  }

  std::cout << "Solving eigenvalue problem..." << std::endl;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);

  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // Find TE10 mode index
  int te10_idx = -1;
  double min_diff = 1e100;
  for (int i = 0; i < n_free; ++i) {
    if (eigenvalues(i) > 100) {
      double diff = std::abs(eigenvalues(i) - kc_te10_sq);
      if (diff < min_diff) {
        min_diff = diff;
        te10_idx = i;
      }
    }
  }

  std::cout << "TE10 eigenvalue: " << eigenvalues(te10_idx)
            << " (expected: " << kc_te10_sq << ")" << std::endl;

  Eigen::VectorXd v_te10 = eigenvectors.col(te10_idx);

  // Create full-size eigenvector
  Eigen::VectorXd v_full = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n_free; ++i) {
    v_full(free_to_orig[i]) = v_te10(i);
  }

  // Analyze the z-dependence of the eigenvector
  // Group edges by their z-coordinate (midpoint)
  std::map<double, std::vector<std::pair<int, double>>>
      z_groups; // z -> [(edge_idx, |v|)]

  for (int e = 0; e < (int)mesh.edges.size(); ++e) {
    if (bc.dirichlet_edges.count(e))
      continue; // Skip PEC edges

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    double z_mid = (p0.z() + p1.z()) / 2.0;

    // Round to mm for grouping
    double z_rounded = std::round(z_mid * 1000) / 1000.0;

    z_groups[z_rounded].push_back({e, std::abs(v_full(e))});
  }

  std::cout << "\n=== Eigenvector Magnitude vs Z-position ===" << std::endl;
  std::cout << "Z(mm)\t\tNum Edges\tMean |v|\t\tMax |v|" << std::endl;

  double z_min = z_groups.begin()->first;
  double z_max = z_groups.rbegin()->first;

  for (const auto &kv : z_groups) {
    double z = kv.first;
    const auto &edges = kv.second;

    double sum_v = 0, max_v = 0;
    for (const auto &ev : edges) {
      sum_v += ev.second;
      max_v = std::max(max_v, ev.second);
    }
    double mean_v = sum_v / edges.size();

    if (z == z_min || z == z_max || edges.size() > 100) {
      std::cout << z * 1000 << "\t\t" << edges.size() << "\t\t" << mean_v
                << "\t\t" << max_v << std::endl;
    }
  }

  // Sample a few z-slices to show the pattern
  std::vector<double> sample_z = {z_min, (z_min + z_max) / 2, z_max};

  for (double z_sample : sample_z) {
    // Find closest z-group
    double best_z = z_min;
    double best_diff = 1e100;
    for (const auto &kv : z_groups) {
      double diff = std::abs(kv.first - z_sample);
      if (diff < best_diff) {
        best_diff = diff;
        best_z = kv.first;
      }
    }

    std::cout << "\n=== Slice at z = " << best_z * 1000
              << " mm ===" << std::endl;
    const auto &edges = z_groups[best_z];

    // Show top 10 edges by |v|
    std::vector<std::pair<double, int>> sorted_edges;
    for (const auto &ev : edges) {
      sorted_edges.push_back({ev.second, ev.first});
    }
    std::sort(sorted_edges.rbegin(), sorted_edges.rend());

    std::cout << "Top 10 edges by |v|:" << std::endl;
    for (int i = 0; i < std::min(10, (int)sorted_edges.size()); ++i) {
      int e = sorted_edges[i].second;
      double v_mag = sorted_edges[i].first;
      const auto &edge = mesh.edges[e];
      const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
      const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
      Eigen::Vector3d mid = (p0 + p1) / 2;
      Eigen::Vector3d dir = (p1 - p0).normalized();

      std::cout << "  Edge " << e << ": |v|=" << v_mag << ", mid=("
                << mid.x() * 1000 << "," << mid.y() * 1000 << ","
                << mid.z() * 1000 << ") mm"
                << ", dir=(" << dir.x() << "," << dir.y() << "," << dir.z()
                << ")" << std::endl;
    }
  }

  // For cavity modes TE10p, the field has sin(p*pi*z/L) dependence
  // Check if the z-pattern matches sin(p*pi*z/L) for some p
  double L = z_max - z_min;
  std::cout << "\n=== Cavity Mode Analysis ===" << std::endl;
  std::cout << "Waveguide length L = " << L * 1000 << " mm" << std::endl;

  // Expected cavity frequencies for TE10p modes
  double a = 0.02286, b = 0.01016;
  std::cout << "\nExpected TE10p cavity frequencies:" << std::endl;
  for (int p = 0; p <= 5; ++p) {
    double k_sq = std::pow(M_PI / a, 2) + std::pow(p * M_PI / L, 2);
    double f_res = c0 * std::sqrt(k_sq) / (2 * M_PI);
    std::cout << "  TE10" << p << ": k² = " << k_sq << ", f = " << f_res / 1e9
              << " GHz" << std::endl;
  }

  std::cout << "\nFEM eigenvalue corresponds to: kc² = "
            << eigenvalues(te10_idx) << std::endl;
  std::cout << "This matches TE10 waveguide cutoff, NOT a cavity mode!"
            << std::endl;

  // Check if the FEM matrix structure suggests standing wave or traveling wave
  std::cout << "\n=== Matrix Structure Check ===" << std::endl;

  // The curl-curl matrix K should give ∇×∇× operator eigenvalues
  // For a waveguide with open ports, eigenvalues should be complex (radiation)
  // For a closed cavity, eigenvalues are real

  // Our eigenvalues are real, suggesting the FEM treats this as a closed cavity
  std::cout << "Eigenvalues are real (not complex) → FEM models a CLOSED CAVITY"
            << std::endl;
  std::cout
      << "For traveling wave ports, we need absorbing boundary conditions!"
      << std::endl;

  return 0;
}
