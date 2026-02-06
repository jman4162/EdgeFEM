// Debug script to investigate FEM core behavior
// Focus on understanding why K dominates M at RF frequencies

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "edgefem/edge_basis.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

namespace {
constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);
} // namespace

int main() {
  std::cout << std::setprecision(8);

  // Physical constants check
  std::cout << "=== Physical Constants ===" << std::endl;
  std::cout << "c0 = " << c0 << " m/s" << std::endl;
  std::cout << "mu0 = " << mu0 << " H/m" << std::endl;
  std::cout << "eps0 = " << eps0 << " F/m" << std::endl;
  std::cout << "mu0 * c0 = " << mu0 * c0 << " (eta0, should be ~377)"
            << std::endl;
  std::cout << "1/(mu0*eps0) = " << 1.0 / (mu0 * eps0) << " (should be c0²)"
            << std::endl;
  std::cout << "c0² = " << c0 * c0 << std::endl;
  std::cout << std::endl;

  // At 10 GHz
  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double lambda = c0 / freq;

  std::cout << "=== At 10 GHz ===" << std::endl;
  std::cout << "omega = " << omega << " rad/s" << std::endl;
  std::cout << "k0 = " << k0 << " rad/m" << std::endl;
  std::cout << "k0² = " << k0 * k0 << std::endl;
  std::cout << "lambda = " << lambda * 1000 << " mm" << std::endl;
  std::cout << std::endl;

  // Element analysis for different sizes
  std::cout << "=== Element Size Analysis ===" << std::endl;
  std::vector<double> sizes = {0.003, 0.001, 0.0003, 0.0001};

  for (double h : sizes) {
    std::array<Eigen::Vector3d, 4> tet = {
        Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(h, 0.0, 0.0),
        Eigen::Vector3d(0.0, h, 0.0), Eigen::Vector3d(0.0, 0.0, h)};

    double V = h * h * h / 6.0; // Volume of regular tet

    auto K = whitney_curl_curl_matrix(tet);
    auto M = whitney_mass_matrix(tet);

    double K_trace = K.trace();
    double M_trace = M.trace();

    // For Whitney elements on a regular tet:
    // K scales as 1/h (curl = O(1/h), integral O(h³), so curl² integral = O(h))
    // M scales as h³ (basis = O(1), integral O(h³))
    // Therefore K/M scales as 1/h⁴

    double K_over_M = K_trace / M_trace;
    double k0_sq = k0 * k0;

    std::cout << "h = " << h * 1000 << " mm:" << std::endl;
    std::cout << "  V = " << V << " m³" << std::endl;
    std::cout << "  K_trace = " << K_trace << ", M_trace = " << M_trace
              << std::endl;
    std::cout << "  K/M = " << K_over_M << std::endl;
    std::cout << "  k0²*M_trace = " << k0_sq * M_trace << std::endl;
    std::cout << "  K/(k0²*M) = " << K_trace / (k0_sq * M_trace) << std::endl;
    std::cout << "  Elements per wavelength: " << lambda / h << std::endl;

    // Generalized eigenvalue problem: K x = λ M x
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K, M);
    double lambda_min = 0;
    for (int i = 0; i < 6; ++i) {
      if (ges.eigenvalues()(i) > 1e-6) {
        lambda_min = ges.eigenvalues()(i);
        break;
      }
    }
    double f_min = std::sqrt(lambda_min) * c0 / (2 * M_PI);
    std::cout << "  Min element eigenfreq: " << f_min / 1e9 << " GHz"
              << std::endl;
    std::cout << std::endl;
  }

  // The key insight: K/M ratio determines element eigenfrequency
  // For an element to resolve frequency f, we need k0² ≈ K/M
  // This means h ∝ 1/√f, not h ∝ 1/f
  //
  // Actually, let's check the scaling more carefully:
  // For regular tet with edge h:
  // - Gradients of barycentric coords scale as 1/h
  // - Curl of Whitney basis = 2 * (grad_a × grad_b) scales as 1/h²
  // - K_ij = V * curl_i · curl_j scales as h³ * (1/h²)² = 1/h
  // - M_ij involves integrals of W_i · W_j
  // - Whitney W = λ_a * grad λ_b - λ_b * grad λ_a
  // - Each term has grad scaling as 1/h, λ as unitless
  // - So W scales as 1/h, and M_ij = ∫ W_i · W_j dV scales as h³ * (1/h)² = h
  //
  // Therefore: K ~ 1/h, M ~ h, and K/M ~ 1/h²
  // Element eigenfrequency: f ~ √(K/M) ~ 1/h
  // So f_element ∝ 1/h, which means smaller elements resolve lower frequencies

  std::cout << "=== Scaling Analysis ===" << std::endl;
  double h1 = 0.003, h2 = 0.001;
  std::array<Eigen::Vector3d, 4> tet1 = {
      Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(h1, 0, 0),
      Eigen::Vector3d(0, h1, 0), Eigen::Vector3d(0, 0, h1)};
  std::array<Eigen::Vector3d, 4> tet2 = {
      Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(h2, 0, 0),
      Eigen::Vector3d(0, h2, 0), Eigen::Vector3d(0, 0, h2)};

  auto K1 = whitney_curl_curl_matrix(tet1);
  auto K2 = whitney_curl_curl_matrix(tet2);
  auto M1 = whitney_mass_matrix(tet1);
  auto M2 = whitney_mass_matrix(tet2);

  std::cout << "h1/h2 = " << h1 / h2 << std::endl;
  std::cout << "K1(0,0)/K2(0,0) = " << K1(0, 0) / K2(0, 0)
            << " (expected h2/h1 = " << h2 / h1 << ")" << std::endl;
  std::cout << "M1(0,0)/M2(0,0) = " << M1(0, 0) / M2(0, 0)
            << " (expected h1/h2 = " << h1 / h2 << ")" << std::endl;

  // So K scales as 1/h (confirmed: 3mm/1mm = 3, K ratio = 0.333)
  // And M scales as h (confirmed: 3mm/1mm = 3, M ratio = 3)

  std::cout << "\n=== Required Element Size ===" << std::endl;
  // For 10 GHz operation, we need element eigenfrequency < 10 GHz
  // Element eigenfrequency ~ 100 GHz for 3mm element
  // Since f_elem ~ 1/h, to get f_elem = 10 GHz from 100 GHz, need h * 10
  // So need 30mm elements? That doesn't make sense physically...

  // Wait - the issue is different. The element eigenfrequency tells us the
  // highest frequency the element can accurately represent, not the lowest.
  // Actually, for wave propagation, we need the wavelength to be resolved
  // by multiple elements, not the frequency to match element eigenfreq.

  // Let me re-examine the actual problem...

  // Load the mesh and check actual assembled matrix eigenvalues
  std::cout << "\n=== Full Mesh Analysis ===" << std::endl;
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  std::cout << "Mesh: " << mesh.nodes.size() << " nodes, " << mesh.edges.size()
            << " edges, " << mesh.tets.size() << " tets" << std::endl;
  std::cout << "PEC edges: " << bc.dirichlet_edges.size() << std::endl;
  std::cout << "Free edges: " << mesh.edges.size() - bc.dirichlet_edges.size()
            << std::endl;

  // Check mesh element sizes
  double min_vol = 1e10, max_vol = 0;
  double total_vol = 0;
  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }
    Eigen::Matrix3d B;
    B.col(0) = X[0] - X[3];
    B.col(1) = X[1] - X[3];
    B.col(2) = X[2] - X[3];
    double V = std::abs(B.determinant()) / 6.0;
    min_vol = std::min(min_vol, V);
    max_vol = std::max(max_vol, V);
    total_vol += V;
  }
  double avg_vol = total_vol / mesh.tets.size();
  double avg_h = std::cbrt(6.0 * avg_vol); // Characteristic length

  std::cout << "Element volumes: min=" << min_vol << ", max=" << max_vol
            << ", avg=" << avg_vol << std::endl;
  std::cout << "Characteristic element size h ≈ " << avg_h * 1000 << " mm"
            << std::endl;
  std::cout << "Elements per wavelength (at 10 GHz): " << lambda / avg_h
            << std::endl;

  // The real issue: Let's check the port formulation
  std::cout << "\n=== Port Weight Analysis ===" << std::endl;

  // Check WR-90 port dimensions
  double WG_A = 0.02286, WG_B = 0.01016;
  double fc = c0 / (2 * WG_A);
  double beta = k0 * std::sqrt(1.0 - std::pow(fc / freq, 2));
  double Z0 = (mu0 * c0) / std::sqrt(1.0 - std::pow(fc / freq, 2));

  std::cout << "WR-90: a=" << WG_A * 1000 << "mm, b=" << WG_B * 1000 << "mm"
            << std::endl;
  std::cout << "fc = " << fc / 1e9 << " GHz" << std::endl;
  std::cout << "beta = " << beta << " rad/m" << std::endl;
  std::cout << "Z0 = " << Z0 << " ohm" << std::endl;

  // Port area
  double port_area = WG_A * WG_B;
  std::cout << "Port area = " << port_area * 1e6 << " mm²" << std::endl;

  // The port admittance added to the matrix is Y = ww^H/Z0
  // For proper matching, we need ||w||² / Z0 to be comparable to matrix
  // diagonal

  // Build ports
  RectWaveguidePort dims{WG_A, WG_B};
  PortSurfaceMesh port1_surface = extract_surface_mesh(mesh, 2);
  PortMode mode = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surface, dims, mode);
  WavePort wp1 = build_wave_port(mesh, port1_surface, mode);

  double w_norm_sq = 0;
  for (int i = 0; i < wp1.weights.size(); ++i) {
    w_norm_sq += std::norm(wp1.weights[i]);
  }
  std::cout << "\nPort weight ||w||² = " << w_norm_sq << std::endl;
  std::cout << "Y_port = ||w||²/Z0 = " << w_norm_sq / Z0 << std::endl;

  // Compare to matrix diagonal
  MaxwellParams p;
  p.omega = omega;
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);

  double avg_diag = 0;
  int count = 0;
  for (int i : wp1.edges) {
    if (!bc.dirichlet_edges.count(i)) {
      avg_diag += std::real(asmbl.A.coeff(i, i));
      count++;
    }
  }
  avg_diag /= count;

  std::cout << "Avg matrix diagonal at port edges: " << avg_diag << std::endl;
  std::cout << "Y_port / A_diag ratio: " << (w_norm_sq / Z0) / avg_diag
            << std::endl;

  // For proper port matching, Y_port should be comparable to A_diag
  // If Y_port << A_diag, port has no effect
  // If Y_port >> A_diag, port dominates

  // The source term is b = 2*w/sqrt(Z0)
  double source_mag = 2.0 / std::sqrt(Z0);
  std::cout << "\nSource scaling 2/sqrt(Z0) = " << source_mag << std::endl;
  std::cout << "Source magnitude (first w): "
            << source_mag * std::abs(wp1.weights[0]) << std::endl;

  // Expected solution scale for matched port:
  // If A*x = b with Y contribution, and V_inc = sqrt(Z0)
  // Then V = w^H * x should be ~ V_inc for matched case
  // This means x ~ V_inc / ||w|| = sqrt(Z0) / ||w||
  double expected_x = std::sqrt(Z0) / std::sqrt(w_norm_sq);
  std::cout << "Expected |x| for matching: " << expected_x << std::endl;

  // The fundamental issue might be that the FEM matrix scale is wrong
  // Let's check what the matrix entries actually represent physically

  std::cout << "\n=== Matrix Entry Physical Interpretation ===" << std::endl;

  // The FEM equation is: (1/μr) ∇×∇×E - k0² εr E = -jωJ
  // In matrix form: (K/μr - k0²εr M) x = b
  // where x are edge DOFs (line integrals of E)
  // and K, M come from element matrices

  // For μr = εr = 1:
  // A = K - k0² M

  // Sample a non-PEC interior edge
  int sample_edge = -1;
  for (int i = 0; i < (int)mesh.edges.size(); ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      // Check if this edge is interior (not on any surface)
      bool on_surface = false;
      for (const auto &tri : mesh.tris) {
        for (int e = 0; e < 3; ++e) {
          if (tri.edges[e] == i) {
            on_surface = true;
            break;
          }
        }
        if (on_surface)
          break;
      }
      if (!on_surface) {
        sample_edge = i;
        break;
      }
    }
  }

  if (sample_edge >= 0) {
    std::complex<double> A_ii = asmbl.A.coeff(sample_edge, sample_edge);
    std::cout << "Sample interior edge " << sample_edge << ": A_ii = " << A_ii
              << std::endl;

    // Estimate K_ii and M_ii contributions
    // A_low = K, A_high = K - k0_high² M
    MaxwellParams p_low;
    p_low.omega = 1e-10;
    auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);
    double K_ii = std::real(asmbl_low.A.coeff(sample_edge, sample_edge));
    double M_ii = (K_ii - std::real(A_ii)) / (k0 * k0);

    std::cout << "Estimated K_ii = " << K_ii << std::endl;
    std::cout << "Estimated M_ii = " << M_ii << std::endl;
    std::cout << "k0² * M_ii = " << k0 * k0 * M_ii << std::endl;
    std::cout << "K_ii / (k0² * M_ii) = " << K_ii / (k0 * k0 * M_ii)
              << std::endl;
  }

  return 0;
}
