// Compare 3D FEM eigenvector pattern at port edges vs expected 2D pattern
// Key question: WHY does the 3D eigenvector have only 3% correlation with 2D mode?

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
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

  double kc_te10_sq = std::pow(M_PI / 0.02286, 2);

  // === 3D FEM Eigensolve ===
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
      double K_ij = std::real(asmbl_low.A.coeff(oi, oj));
      double A_test_ij = std::real(asmbl_test.A.coeff(oi, oj));
      K_free(i, j) = K_ij;
      M_free(i, j) = (K_ij - A_test_ij) / (k0_test * k0_test);
    }
  }

  std::cout << "Solving 3D eigenvalue problem (" << n_free << " DOFs)..." << std::endl;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);
  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // Find TE10
  int te10_idx = -1;
  double min_diff = 1e100;
  for (int i = 0; i < n_free; ++i) {
    if (eigenvalues(i) > 100) {
      double diff = std::abs(eigenvalues(i) - kc_te10_sq);
      if (diff < min_diff) { min_diff = diff; te10_idx = i; }
    }
  }

  std::cout << "TE10 eigenvalue: " << eigenvalues(te10_idx) << " (expected: " << kc_te10_sq << ")" << std::endl;

  // Create full eigenvector
  Eigen::VectorXd v_3d = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n_free; ++i) {
    v_3d(free_to_orig[i]) = eigenvectors.col(te10_idx)(i);
  }

  // === Get port 1 surface mesh and edges ===
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  RectWaveguidePort dims{0.02286, 0.01016};
  PortMode mode1 = solve_te10_mode(dims, 10e9);
  populate_te10_field(port1_surf, dims, mode1);
  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);

  std::cout << "\n=== Port 1 Edge Analysis ===" << std::endl;
  std::cout << "Total port edges: " << wp1.edges.size() << std::endl;

  // For each free port edge, compare:
  // 1. Analytical weight w_e (from mode1)
  // 2. 3D FEM eigenvector v_e
  // 3. Expected pattern based on edge position

  std::cout << "\nEdge analysis (sorted by |w|):" << std::endl;
  std::cout << "Edge\t\tMidpoint(mm)\t\t|w|\t\t|v_3d|\t\tratio" << std::endl;

  std::vector<std::tuple<int, double, double, double, Eigen::Vector3d>> edge_data;

  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    Eigen::Vector3d mid = (p0 + p1) / 2;

    double w_mag = std::abs(wp1.weights(i));
    double v_mag = std::abs(v_3d(e));

    edge_data.push_back({e, w_mag, v_mag, (w_mag > 1e-10 ? v_mag/w_mag : 0), mid});
  }

  // Sort by |w| descending
  std::sort(edge_data.begin(), edge_data.end(),
            [](const auto &a, const auto &b) { return std::get<1>(a) > std::get<1>(b); });

  // Show top 20
  for (int i = 0; i < std::min(20, (int)edge_data.size()); ++i) {
    auto &[e, w_mag, v_mag, ratio, mid] = edge_data[i];
    std::cout << e << "\t\t(" << mid.x()*1000 << "," << mid.y()*1000 << ")\t"
              << w_mag << "\t" << v_mag << "\t" << ratio << std::endl;
  }

  // Check if the edge direction matters
  std::cout << "\n=== Edge Direction Analysis ===" << std::endl;
  std::cout << "TE10 has Ey only, so y-directed edges should dominate" << std::endl;

  double sum_y_edges_w = 0, sum_x_edges_w = 0;
  double sum_y_edges_v = 0, sum_x_edges_v = 0;
  int n_y_edges = 0, n_x_edges = 0;

  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    Eigen::Vector3d dir = (p1 - p0).normalized();

    double w_mag = std::abs(wp1.weights(i));
    double v_mag = std::abs(v_3d(e));

    // Classify edge as x-directed or y-directed
    if (std::abs(dir.x()) > std::abs(dir.y())) {
      sum_x_edges_w += w_mag;
      sum_x_edges_v += v_mag;
      n_x_edges++;
    } else {
      sum_y_edges_w += w_mag;
      sum_y_edges_v += v_mag;
      n_y_edges++;
    }
  }

  std::cout << "X-directed edges: " << n_x_edges << ", sum |w| = " << sum_x_edges_w << ", sum |v_3d| = " << sum_x_edges_v << std::endl;
  std::cout << "Y-directed edges: " << n_y_edges << ", sum |w| = " << sum_y_edges_w << ", sum |v_3d| = " << sum_y_edges_v << std::endl;
  std::cout << "Ratio (Y/X) for w: " << sum_y_edges_w / sum_x_edges_w << std::endl;
  std::cout << "Ratio (Y/X) for v_3d: " << sum_y_edges_v / sum_x_edges_v << std::endl;

  // Check if the 3D eigenvector has significant z-variation at port edges
  std::cout << "\n=== 3D Eigenvector Pattern ===" << std::endl;

  // For true waveguide mode, v_3d should be constant along z at a given (x,y)
  // For cavity mode, v_3d has sin(pπz/L) variation

  // Group edges by approximate (x,y) position and see if z affects value
  std::map<std::pair<int,int>, std::vector<std::pair<double, double>>> xy_to_z_v;

  for (int e = 0; e < (int)mesh.edges.size(); ++e) {
    if (bc.dirichlet_edges.count(e)) continue;

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    Eigen::Vector3d mid = (p0 + p1) / 2;

    // Round to mm for grouping
    int x_mm = (int)(mid.x() * 100);
    int y_mm = (int)(mid.y() * 100);

    xy_to_z_v[{x_mm, y_mm}].push_back({mid.z(), std::abs(v_3d(e))});
  }

  // Find an (x,y) with many z values
  int max_z_count = 0;
  std::pair<int,int> best_xy;
  for (const auto &kv : xy_to_z_v) {
    if ((int)kv.second.size() > max_z_count) {
      max_z_count = kv.second.size();
      best_xy = kv.first;
    }
  }

  if (max_z_count > 3) {
    std::cout << "Eigenvector at (x,y)≈(" << best_xy.first/100.0 << "," << best_xy.second/100.0 << ") mm:" << std::endl;
    auto &z_vals = xy_to_z_v[best_xy];
    std::sort(z_vals.begin(), z_vals.end());
    std::cout << "z(mm)\t\t|v_3d|" << std::endl;
    for (const auto &zv : z_vals) {
      std::cout << zv.first*1000 << "\t\t" << zv.second << std::endl;
    }
  }

  // Finally, check at port edges specifically
  std::cout << "\n=== At Port 1 (z≈0) Specifically ===" << std::endl;

  // For TE10: Ey ∝ sin(πx/a), which varies from 0 at x=0 to 1 at x=a/2 to 0 at x=a
  // So edges near x=a/2 ≈ 11.4mm should have largest weight

  for (int i = 0; i < std::min(10, (int)edge_data.size()); ++i) {
    auto &[e, w_mag, v_mag, ratio, mid] = edge_data[i];
    double x = mid.x();
    double expected_sin = std::sin(M_PI * x / 0.02286);
    std::cout << "Edge " << e << " at x=" << x*1000 << "mm: "
              << "|w|=" << w_mag << " (expect ~sin=" << expected_sin << "), "
              << "|v_3d|=" << v_mag << std::endl;
  }

  return 0;
}
