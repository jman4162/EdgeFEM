// Test solving 2D eigenvalue problem on port surface
// This gives the correct discrete mode pattern for port weights

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;

// Build 2D curl-curl and mass matrices for port surface
// For TE modes, we solve: ∇t²Hz + kc²Hz = 0 (scalar Helmholtz on surface)
// With Neumann BC: ∂Hz/∂n = 0 on PEC walls
void build_2d_fem_matrices(const Mesh &mesh, const PortSurfaceMesh &surface,
                           const BC &bc, Eigen::MatrixXd &K,
                           Eigen::MatrixXd &M) {
  // Get node count for surface mesh
  int num_nodes = surface.mesh.nodes.size();

  K = Eigen::MatrixXd::Zero(num_nodes, num_nodes);
  M = Eigen::MatrixXd::Zero(num_nodes, num_nodes);

  // Assemble over surface triangles
  for (size_t tri_idx = 0; tri_idx < surface.mesh.tris.size(); ++tri_idx) {
    const auto &tri2d = surface.mesh.tris[tri_idx];

    // Get node positions (using 2D coordinates from surface.mesh)
    Eigen::Vector2d p0(
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[0])].xyz.x(),
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[0])].xyz.y());
    Eigen::Vector2d p1(
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[1])].xyz.x(),
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[1])].xyz.y());
    Eigen::Vector2d p2(
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[2])].xyz.x(),
        surface.mesh.nodes[surface.mesh.nodeIndex.at(tri2d.conn[2])].xyz.y());

    // Triangle area
    double area = 0.5 * std::abs((p1.x() - p0.x()) * (p2.y() - p0.y()) -
                                 (p2.x() - p0.x()) * (p1.y() - p0.y()));

    if (area < 1e-15)
      continue;

    // Shape function gradients for linear triangular elements
    // For node i in triangle: N_i = (a_i + b_i*x + c_i*y) / (2*area)
    // grad N_i = (b_i, c_i) / (2*area)
    Eigen::Matrix<double, 3, 2>
        grad_N; // Each row is gradient of one shape function

    double b0 = p1.y() - p2.y();
    double c0 = p2.x() - p1.x();
    double b1 = p2.y() - p0.y();
    double c1 = p0.x() - p2.x();
    double b2 = p0.y() - p1.y();
    double c2 = p1.x() - p0.x();

    grad_N(0, 0) = b0 / (2 * area);
    grad_N(0, 1) = c0 / (2 * area);
    grad_N(1, 0) = b1 / (2 * area);
    grad_N(1, 1) = c1 / (2 * area);
    grad_N(2, 0) = b2 / (2 * area);
    grad_N(2, 1) = c2 / (2 * area);

    // Local stiffness: K_ij = ∫ (∇N_i · ∇N_j) dA = area * (∇N_i · ∇N_j)
    // Local mass: M_ij = ∫ N_i N_j dA = area/12 * (1 + δ_ij)
    for (int i = 0; i < 3; ++i) {
      int gi = surface.mesh.nodeIndex.at(tri2d.conn[i]);
      for (int j = 0; j < 3; ++j) {
        int gj = surface.mesh.nodeIndex.at(tri2d.conn[j]);

        // Stiffness contribution
        double k_ij = area * grad_N.row(i).dot(grad_N.row(j));
        K(gi, gj) += k_ij;

        // Mass contribution (lumped mass would be area/3 on diagonal)
        double m_ij = area * (i == j ? 1.0 / 6.0 : 1.0 / 12.0);
        M(gi, gj) += m_ij;
      }
    }
  }
}

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2 * M_PI * freq;

  RectWaveguidePort dims{0.02286, 0.01016};

  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);

  std::cout << "=== 2D Port Surface Mesh ===" << std::endl;
  std::cout << "Nodes: " << port1_surf.mesh.nodes.size() << std::endl;
  std::cout << "Triangles: " << port1_surf.mesh.tris.size() << std::endl;
  std::cout << "Boundary lines: " << port1_surf.mesh.boundary_lines.size()
            << std::endl;

  // Build 2D FEM matrices
  Eigen::MatrixXd K_2d, M_2d;
  build_2d_fem_matrices(mesh, port1_surf, bc, K_2d, M_2d);

  int n = K_2d.rows();
  std::cout << "\n2D FEM matrices: " << n << " x " << n << std::endl;
  std::cout << "K trace = " << K_2d.trace() << std::endl;
  std::cout << "M trace = " << M_2d.trace() << std::endl;

  // Solve 2D eigenvalue problem: K*Hz = kc² * M * Hz
  // For Neumann BC on all boundaries, we don't eliminate any DOFs
  std::cout << "\nSolving 2D eigenvalue problem..." << std::endl;

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_2d, M_2d);

  if (ges.info() != Eigen::Success) {
    std::cout << "Eigenvalue solver failed!" << std::endl;
    return 1;
  }

  auto eigenvalues = ges.eigenvalues();
  auto eigenvectors = ges.eigenvectors();

  // Expected kc² values for rectangular waveguide TE modes
  double kc_te10_sq = std::pow(M_PI / dims.a, 2);     // TE10
  double kc_te20_sq = std::pow(2 * M_PI / dims.a, 2); // TE20
  double kc_te01_sq = std::pow(M_PI / dims.b, 2);     // TE01

  std::cout << "\nExpected kc² values:" << std::endl;
  std::cout << "  TE10: " << kc_te10_sq << std::endl;
  std::cout << "  TE20: " << kc_te20_sq << std::endl;
  std::cout << "  TE01: " << kc_te01_sq << std::endl;

  std::cout << "\nFirst 10 eigenvalues (kc²):" << std::endl;
  for (int i = 0; i < std::min(10, n); ++i) {
    std::cout << "  λ_" << i << " = " << eigenvalues(i);
    if (std::abs(eigenvalues(i) - kc_te10_sq) / kc_te10_sq < 0.1)
      std::cout << " ← TE10";
    else if (std::abs(eigenvalues(i) - kc_te20_sq) / kc_te20_sq < 0.1)
      std::cout << " ← TE20";
    else if (std::abs(eigenvalues(i) - kc_te01_sq) / kc_te01_sq < 0.1)
      std::cout << " ← TE01";
    else if (eigenvalues(i) < 1.0)
      std::cout << " ← DC mode";
    std::cout << std::endl;
  }

  // Find TE10 mode index
  int te10_idx = -1;
  double min_diff = 1e100;
  for (int i = 0; i < n; ++i) {
    double diff = std::abs(eigenvalues(i) - kc_te10_sq);
    if (diff < min_diff) {
      min_diff = diff;
      te10_idx = i;
    }
  }

  std::cout << "\nTE10 mode: index " << te10_idx
            << ", kc² = " << eigenvalues(te10_idx)
            << " (error: " << min_diff / kc_te10_sq * 100 << "%)" << std::endl;

  // Get TE10 eigenvector (Hz pattern on port surface)
  Eigen::VectorXd Hz_te10 = eigenvectors.col(te10_idx);

  // Normalize for visualization
  Hz_te10 /= Hz_te10.maxCoeff();

  std::cout << "\nHz pattern at select nodes:" << std::endl;
  // Sample some nodes to verify cos(πx/a) pattern
  for (int i = 0; i < std::min(10, n); ++i) {
    const auto &node = port1_surf.mesh.nodes[i];
    double x = node.xyz.x();
    double expected = std::cos(M_PI * x / dims.a);
    std::cout << "  Node " << i << " at x=" << x * 1000
              << "mm: Hz=" << Hz_te10(i) << ", expected cos=" << expected
              << std::endl;
  }

  // Now compute port weights using 2D FEM eigenvector
  // Et = jωμ/kc² × (ẑ × ∇t Hz) = jωμ/kc² × (∂Hz/∂y, -∂Hz/∂x, 0)
  // For TE10 with Hz = cos(πx/a): ∂Hz/∂x = -π/a sin(πx/a), ∂Hz/∂y = 0
  // So Et = jωμ/kc² × (0, π/a sin(πx/a), 0) = jωμ/(akc²) × (0, π sin(πx/a), 0)

  // For the FEM eigenvector, we compute ∇Hz at each triangle, then integrate Et
  // along edges
  std::cout << "\n=== Computing FEM-based Port Weights ===" << std::endl;

  PortMode mode1 = solve_te10_mode(dims, freq);
  double kc = mode1.kc;
  double mu = std::real(mode1.mu);
  std::complex<double> j_factor(0, omega * mu / (kc * kc));

  // Create weight vector for port edges using analytical Hz
  populate_te10_field(port1_surf, dims, mode1);
  WavePort wp_analytical = build_wave_port(mesh, port1_surf, mode1);

  // Now create weights using FEM eigenvector Hz
  // The FEM eigenvector has entries for nodes in surface.mesh (2D indices)
  // We need to populate mode1_fem.field with values at the correct indices

  PortMode mode1_fem = mode1;
  mode1_fem.field.resize(port1_surf.mesh.nodes.size());

  // Use the FEM kc from eigenvalue
  mode1_fem.kc = std::sqrt(eigenvalues(te10_idx));

  // The surface mesh node indices go from 0 to n-1
  for (int i = 0; i < n; ++i) {
    // Scale the FEM eigenvector to match expected amplitude
    // Analytical Hz uses unit-power normalization with some amplitude A
    // FEM Hz_te10 is an eigenvector (normalized to unit norm typically)
    // We'll scale FEM Hz to have similar range as analytical Hz
    mode1_fem.field(i) = Hz_te10(i) * mode1.field.cwiseAbs().maxCoeff();
  }

  // Build wave port with FEM Hz field
  WavePort wp_fem_actual = build_wave_port(mesh, port1_surf, mode1_fem);

  std::cout << "Analytical weights ||w||² = "
            << wp_analytical.weights.squaredNorm() << std::endl;
  std::cout << "FEM-based weights ||w||² = "
            << wp_fem_actual.weights.squaredNorm() << std::endl;

  // Check correlation between analytical and FEM weights
  Eigen::VectorXcd w_ana = wp_analytical.weights;
  Eigen::VectorXcd w_fem = wp_fem_actual.weights;

  std::complex<double> corr = (w_ana.conjugate().transpose() * w_fem)(0, 0);
  corr /= (w_ana.norm() * w_fem.norm());

  std::cout << "Correlation |<w_ana, w_fem>| = " << std::abs(corr) << std::endl;

  // Compare weight patterns
  std::cout << "\nFirst 10 weight comparisons:" << std::endl;
  std::cout << "Edge\t\tw_analytical\t\tw_fem\t\t\tratio" << std::endl;
  int shown = 0;
  for (size_t i = 0; i < wp_analytical.edges.size() && shown < 10; ++i) {
    int e = wp_analytical.edges[i];
    if (!bc.dirichlet_edges.count(e)) {
      auto w_a = wp_analytical.weights(i);
      auto w_f = wp_fem_actual.weights(i);
      auto ratio =
          (std::abs(w_a) > 1e-10) ? w_f / w_a : std::complex<double>(0);
      std::cout << e << "\t\t" << w_a << "\t\t" << w_f << "\t\t" << ratio
                << std::endl;
      shown++;
    }
  }

  // Now test S-parameters with FEM-based weights
  std::cout << "\n=== S-parameters Test ===" << std::endl;

  // Normalize FEM weights for proper S-param extraction
  double Z0 = std::real(mode1.Z0);
  double target_norm_sq = std::sqrt(Z0);

  double fem_norm_sq = wp_fem_actual.weights.squaredNorm();
  double scale_fem = std::sqrt(target_norm_sq / fem_norm_sq);

  WavePort wp1_fem_norm = wp_fem_actual;
  wp1_fem_norm.weights *= scale_fem;
  wp1_fem_norm.mode = mode1;

  // Also create port 2 with FEM weights
  // We need to solve 2D eigenvalue on port 2 surface as well
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Build 2D matrices for port 2
  Eigen::MatrixXd K_2d_p2, M_2d_p2;
  build_2d_fem_matrices(mesh, port2_surf, bc, K_2d_p2, M_2d_p2);

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges_p2(K_2d_p2,
                                                                   M_2d_p2);
  auto eigenvalues_p2 = ges_p2.eigenvalues();
  auto eigenvectors_p2 = ges_p2.eigenvectors();

  // Find TE10 mode for port 2
  int te10_idx_p2 = -1;
  double min_diff_p2 = 1e100;
  for (int i = 0; i < (int)eigenvalues_p2.size(); ++i) {
    double diff = std::abs(eigenvalues_p2(i) - kc_te10_sq);
    if (diff < min_diff_p2) {
      min_diff_p2 = diff;
      te10_idx_p2 = i;
    }
  }
  Eigen::VectorXd Hz_te10_p2 = eigenvectors_p2.col(te10_idx_p2);
  Hz_te10_p2 /= Hz_te10_p2.maxCoeff();

  PortMode mode2_fem = mode1;
  mode2_fem.field.resize(port2_surf.mesh.nodes.size());
  mode2_fem.kc = std::sqrt(eigenvalues_p2(te10_idx_p2));

  for (int i = 0; i < (int)port2_surf.mesh.nodes.size(); ++i) {
    mode2_fem.field(i) = Hz_te10_p2(i) * mode1.field.cwiseAbs().maxCoeff();
  }

  WavePort wp2_fem = build_wave_port(mesh, port2_surf, mode2_fem);
  double fem_norm_sq_p2 = wp2_fem.weights.squaredNorm();
  double scale_fem_p2 = std::sqrt(target_norm_sq / fem_norm_sq_p2);
  wp2_fem.weights *= scale_fem_p2;
  wp2_fem.mode = mode1;

  std::cout << "Normalized ||w||² = " << wp1_fem_norm.weights.squaredNorm()
            << " (target: " << target_norm_sq << ")" << std::endl;

  MaxwellParams p;
  p.omega = omega;

  std::vector<WavePort> ports_fem{wp1_fem_norm, wp2_fem};
  auto S_fem = calculate_sparams(mesh, p, bc, ports_fem);

  std::cout << "S11 = " << S_fem(0, 0) << " (|S11| = " << std::abs(S_fem(0, 0))
            << ")" << std::endl;
  std::cout << "S21 = " << S_fem(1, 0) << " (|S21| = " << std::abs(S_fem(1, 0))
            << ")" << std::endl;

  double beta = std::real(mode1.beta);
  double L = 0.05;
  std::cout << "\nExpected: S11 ≈ 0, S21 ≈ exp(-jβL) = "
            << std::exp(std::complex<double>(0, -beta * L)) << std::endl;

  return 0;
}
