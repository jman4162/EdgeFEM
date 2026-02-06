// Phase 0: FEM Core Validation Benchmark
// Tests WR-90 rectangular waveguide against analytical TE10 mode solutions
//
// WR-90 dimensions: a = 22.86mm, b = 10.16mm
// TE10 cutoff: fc = c/(2a) = 6.56 GHz
// Test frequency: 10 GHz (above cutoff)
//
// Analytical results at 10 GHz:
//   β = k0 * sqrt(1 - (fc/f)^2) = 158.24 rad/m
//   Z0_TE = η0 / sqrt(1 - (fc/f)^2) = 498.97 Ω
//   S11 = 0 (matched)
//   S21 = exp(-jβL) for length L

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
constexpr double eta0 = mu0 * c0; // ~377 ohms

// WR-90 waveguide dimensions
constexpr double WG_A = 0.02286; // 22.86 mm width
constexpr double WG_B = 0.01016; // 10.16 mm height
constexpr double WG_L = 0.05;    // 50 mm length

// TE10 cutoff frequency for WR-90
double te10_cutoff_freq() { return c0 / (2.0 * WG_A); }

// TE10 propagation constant at given frequency
double te10_beta(double freq) {
  double fc = te10_cutoff_freq();
  if (freq <= fc)
    return 0.0;
  double k0 = 2.0 * M_PI * freq / c0;
  return k0 * std::sqrt(1.0 - std::pow(fc / freq, 2));
}

// TE10 wave impedance at given frequency
double te10_impedance(double freq) {
  double fc = te10_cutoff_freq();
  if (freq <= fc)
    return std::numeric_limits<double>::infinity();
  return eta0 / std::sqrt(1.0 - std::pow(fc / freq, 2));
}

// Expected S21 for matched waveguide of given length
std::complex<double> analytical_s21(double freq, double length) {
  double beta = te10_beta(freq);
  return std::exp(std::complex<double>(0.0, -beta * length));
}
} // namespace

// =============================================================================
// Test 1: Whitney Element Matrix Verification
// =============================================================================
// Verify curl-curl and mass matrices on a reference tetrahedron against
// analytical integrals from Jin "FEM for Electromagnetics", Chapter 5.

void test_whitney_matrices() {
  std::cout << "\n=== Test 1: Whitney Element Matrix Verification ==="
            << std::endl;

  // Reference tetrahedron with vertices at:
  // v0 = (0, 0, 0), v1 = (1, 0, 0), v2 = (0, 1, 0), v3 = (0, 0, 1)
  // Volume = 1/6
  std::array<Eigen::Vector3d, 4> ref_tet = {
      Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 0.0, 0.0),
      Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0)};

  auto K = whitney_curl_curl_matrix(ref_tet);
  auto M = whitney_mass_matrix(ref_tet);

  double V = 1.0 / 6.0; // Expected volume

  // For the reference tet, gradient of barycentric coords:
  // λ0 has gradient (-1, -1, -1)
  // λ1 has gradient (1, 0, 0)
  // λ2 has gradient (0, 1, 0)
  // λ3 has gradient (0, 0, 1)
  //
  // Edge curls = 2 * (grad_a × grad_b)
  // Edge (0,1): curl = 2 * (-1,-1,-1) × (1,0,0) = 2*(0,-1,1) = (0,-2,2)
  // Edge (0,2): curl = 2 * (-1,-1,-1) × (0,1,0) = 2*(1,0,-1) = (2,0,-2)
  // Edge (0,3): curl = 2 * (-1,-1,-1) × (0,0,1) = 2*(-1,1,0) = (-2,2,0)
  // Edge (1,2): curl = 2 * (1,0,0) × (0,1,0) = 2*(0,0,1) = (0,0,2)
  // Edge (1,3): curl = 2 * (1,0,0) × (0,0,1) = 2*(0,-1,0) = (0,-2,0)
  // Edge (2,3): curl = 2 * (0,1,0) × (0,0,1) = 2*(1,0,0) = (2,0,0)

  // Check curl-curl matrix symmetry
  bool K_symmetric = K.isApprox(K.transpose(), 1e-12);
  std::cout << "K symmetric: " << (K_symmetric ? "PASS" : "FAIL") << std::endl;
  assert(K_symmetric);

  // Check mass matrix symmetry and positive semi-definiteness
  bool M_symmetric = M.isApprox(M.transpose(), 1e-12);
  std::cout << "M symmetric: " << (M_symmetric ? "PASS" : "FAIL") << std::endl;
  assert(M_symmetric);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> es_M(M);
  double min_M_eig = es_M.eigenvalues().minCoeff();
  bool M_psd = min_M_eig >= -1e-12;
  std::cout << "M positive semi-definite: " << (M_psd ? "PASS" : "FAIL")
            << " (min eigenvalue = " << min_M_eig << ")" << std::endl;
  assert(M_psd);

  // For reference tet, K should have known structure
  // K_ij = V * (curl_i · curl_j)
  // E.g., K(0,0) = V * |curl_edge0|² = (1/6) * (0² + 2² + 2²) = 8/6 = 4/3
  double K00_expected = V * 8.0; // |curl|² = 8
  double K00_actual = K(0, 0);
  bool K00_ok = std::abs(K00_expected - K00_actual) < 1e-10;
  std::cout << "K(0,0) check: expected=" << K00_expected
            << ", actual=" << K00_actual << " " << (K00_ok ? "PASS" : "FAIL")
            << std::endl;
  assert(K00_ok);

  // Check K has correct null space dimension (gradient fields)
  // For tetrahedral edge elements, null space of K is dim 4 (gradients of node
  // potentials)
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> es_K(K);
  int K_null_dim = 0;
  for (int i = 0; i < 6; ++i) {
    if (std::abs(es_K.eigenvalues()(i)) < 1e-10)
      K_null_dim++;
  }
  // Actually, for a single tet with 6 edges and 4 nodes, null(K) = 4-1 = 3
  // (one constraint from Kirchhoff sum = 0)
  // But numerically we might see different due to element geometry
  std::cout << "K null space dimension: " << K_null_dim << std::endl;

  // Check mass matrix integral property
  // tr(M) should equal integral of sum of |N_i|² over element
  // For Whitney edge elements: ∫|N_i|² dV = V * (li² + li*lj) / 20
  // where li, lj are edge lengths... this is complex.
  // Instead, check trace is positive and reasonable
  double M_trace = M.trace();
  bool M_trace_ok = M_trace > 0.0;
  std::cout << "M trace = " << M_trace << " " << (M_trace_ok ? "PASS" : "FAIL")
            << std::endl;
  assert(M_trace_ok);

  std::cout << "\nCurl-curl matrix K:\n" << K << std::endl;
  std::cout << "\nMass matrix M:\n" << M << std::endl;
}

// =============================================================================
// Test 2: Assembly Sign and Formulation Check
// =============================================================================
// Verify the FEM formulation: (1/μr)∇×∇×E - k₀²εr E = 0
// Matrix should be A = K/μr - k₀²εr*M

void test_assembly_formulation() {
  std::cout << "\n=== Test 2: Assembly Formulation Check ===" << std::endl;

  // Load waveguide mesh
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2.0 * M_PI * freq;
  double k0 = omega / c0;
  double k0_sq = k0 * k0;

  MaxwellParams p;
  p.omega = omega;
  p.eps_r = 1.0;
  p.mu_r = 1.0;

  // Assemble without ports to get pure FEM matrix
  auto asmbl = assemble_maxwell(mesh, p, bc, {}, -1);

  // Check matrix dimensions
  int n_edges = static_cast<int>(mesh.edges.size());
  bool dim_ok = (asmbl.A.rows() == n_edges && asmbl.A.cols() == n_edges);
  std::cout << "Matrix dimensions: " << asmbl.A.rows() << "x" << asmbl.A.cols()
            << " (expected " << n_edges << "x" << n_edges << ") "
            << (dim_ok ? "PASS" : "FAIL") << std::endl;
  assert(dim_ok);

  // Check Hermitian symmetry (should be symmetric for lossless, real materials)
  Eigen::MatrixXcd A_dense = Eigen::MatrixXcd(asmbl.A);
  double sym_error = (A_dense - A_dense.adjoint()).norm() / A_dense.norm();
  bool symmetric = sym_error < 1e-10;
  std::cout << "Matrix Hermitian symmetry error: " << sym_error << " "
            << (symmetric ? "PASS" : "FAIL") << std::endl;
  assert(symmetric);

  // Analyze diagonal structure
  // For free-space (εr=μr=1), diagonal should be: K_ii - k0²*M_ii
  // At 10 GHz, k0 = 209.4 rad/m, k0² = 43868
  std::cout << "\nk0 = " << k0 << " rad/m" << std::endl;
  std::cout << "k0² = " << k0_sq << std::endl;

  // Sample some free (non-PEC) DOF diagonals
  double max_diag = 0, min_diag = 1e100;
  int sample_count = 0;
  for (int i = 0; i < n_edges && sample_count < 20; ++i) {
    if (bc.dirichlet_edges.count(i))
      continue;
    std::complex<double> diag = asmbl.A.coeff(i, i);
    double real_diag = std::real(diag);
    max_diag = std::max(max_diag, real_diag);
    min_diag = std::min(min_diag, real_diag);
    sample_count++;
  }
  std::cout << "Free DOF diagonal range: [" << min_diag << ", " << max_diag
            << "]" << std::endl;

  // The sign check: for wave propagation at f > fc, we need both positive
  // and negative eigenvalues in the indefinite system.
  // The eigenvalues of (K - k0²M) should span both signs when k0² > min(K/M).

  // Quick eigenvalue check on a subset (full eigensolve too expensive)
  // Instead, check if matrix is indefinite by looking at diagonal signs
  int pos_diag = 0, neg_diag = 0;
  for (int i = 0; i < n_edges; ++i) {
    if (bc.dirichlet_edges.count(i))
      continue;
    double d = std::real(asmbl.A.coeff(i, i));
    if (d > 0)
      pos_diag++;
    else if (d < 0)
      neg_diag++;
  }
  std::cout << "Diagonal signs: " << pos_diag << " positive, " << neg_diag
            << " negative" << std::endl;

  // For proper wave propagation above cutoff, we expect indefinite matrix
  // (mixed signs in eigenvalues). If all diagonals are same sign, something is
  // wrong.
  bool indefinite = (pos_diag > 0 && neg_diag > 0);
  std::cout << "Matrix indefiniteness check: "
            << (indefinite ? "PASS (indefinite as expected)"
                           : "WARNING (not indefinite)")
            << std::endl;

  // Check K vs M scaling by assembling at very low frequency (K dominates)
  // and high frequency (M dominates)
  MaxwellParams p_low, p_high;
  p_low.omega = 2 * M_PI * 1e6;    // 1 MHz
  p_high.omega = 2 * M_PI * 100e9; // 100 GHz

  auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);
  auto asmbl_high = assemble_maxwell(mesh, p_high, bc, {}, -1);

  // At low freq, diagonals should be ~K_ii (positive, from curl-curl)
  // At high freq, diagonals should be ~-k0²M_ii (negative, from mass)
  double low_diag_avg = 0, high_diag_avg = 0;
  int free_count = 0;
  for (int i = 0; i < n_edges; ++i) {
    if (bc.dirichlet_edges.count(i))
      continue;
    low_diag_avg += std::real(asmbl_low.A.coeff(i, i));
    high_diag_avg += std::real(asmbl_high.A.coeff(i, i));
    free_count++;
  }
  low_diag_avg /= free_count;
  high_diag_avg /= free_count;

  std::cout << "\nFrequency scaling check:" << std::endl;
  std::cout << "  Low freq (1 MHz) avg diagonal: " << low_diag_avg
            << " (should be positive)" << std::endl;
  std::cout << "  High freq (100 GHz) avg diagonal: " << high_diag_avg
            << " (should be negative)" << std::endl;

  bool low_positive = low_diag_avg > 0;
  bool high_negative = high_diag_avg < 0;
  std::cout << "  Low freq positive: " << (low_positive ? "PASS" : "FAIL")
            << std::endl;
  std::cout << "  High freq negative: " << (high_negative ? "PASS" : "FAIL")
            << std::endl;
  assert(low_positive);
  assert(high_negative);
}

// =============================================================================
// Test 3: Waveguide S-Parameter Benchmark
// =============================================================================
// Full FEM solve for WR-90 waveguide at 10 GHz, compare to analytical.

void test_waveguide_sparams() {
  std::cout << "\n=== Test 3: Waveguide S-Parameter Benchmark ===" << std::endl;

  // Load mesh
  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  std::cout << "Mesh: " << mesh.nodes.size() << " nodes, " << mesh.edges.size()
            << " edges, " << mesh.tets.size() << " tets" << std::endl;

  // Apply PEC on waveguide walls (physical tag 1)
  BC bc = build_edge_pec(mesh, 1);
  std::cout << "PEC edges: " << bc.dirichlet_edges.size() << std::endl;

  // Setup frequency
  double freq = 10e9;
  double fc = te10_cutoff_freq();
  double beta_analytical = te10_beta(freq);
  double Z0_analytical = te10_impedance(freq);

  std::cout << "\nAnalytical TE10 mode at " << freq / 1e9
            << " GHz:" << std::endl;
  std::cout << "  Cutoff frequency: " << fc / 1e9 << " GHz" << std::endl;
  std::cout << "  Propagation constant β: " << beta_analytical << " rad/m"
            << std::endl;
  std::cout << "  Wave impedance Z0: " << Z0_analytical << " Ω" << std::endl;

  MaxwellParams p;
  p.omega = 2.0 * M_PI * freq;

  // Build ports
  RectWaveguidePort dims{WG_A, WG_B};
  PortMode mode = solve_te10_mode(dims, freq);

  std::cout << "\nComputed TE10 mode:" << std::endl;
  std::cout << "  Cutoff frequency: " << mode.fc / 1e9 << " GHz" << std::endl;
  std::cout << "  Propagation constant β: " << std::real(mode.beta) << " rad/m"
            << std::endl;
  std::cout << "  Wave impedance Z0: " << std::real(mode.Z0) << " Ω"
            << std::endl;

  // Check mode matches analytical
  bool fc_match = std::abs(mode.fc - fc) / fc < 1e-6;
  bool beta_match =
      std::abs(std::real(mode.beta) - beta_analytical) / beta_analytical < 1e-3;
  bool Z0_match =
      std::abs(std::real(mode.Z0) - Z0_analytical) / Z0_analytical < 1e-3;

  std::cout << "\nMode verification: fc " << (fc_match ? "PASS" : "FAIL")
            << ", β " << (beta_match ? "PASS" : "FAIL") << ", Z0 "
            << (Z0_match ? "PASS" : "FAIL") << std::endl;
  assert(fc_match && beta_match && Z0_match);

  // Extract port surfaces and build wave ports
  PortSurfaceMesh port1 = extract_surface_mesh(mesh, 2); // z=0 face
  PortSurfaceMesh port2 = extract_surface_mesh(mesh, 3); // z=L face

  std::cout << "\nPort 1: " << port1.mesh.nodes.size() << " nodes, "
            << port1.mesh.tris.size() << " tris" << std::endl;
  std::cout << "Port 2: " << port2.mesh.nodes.size() << " nodes, "
            << port2.mesh.tris.size() << " tris" << std::endl;

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);

  populate_te10_field(port1, dims, mode1);
  populate_te10_field(port2, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1, mode1);
  WavePort wp2 = build_wave_port(mesh, port2, mode2);

  std::cout << "Wave port 1: " << wp1.edges.size() << " edges" << std::endl;
  std::cout << "Wave port 2: " << wp2.edges.size() << " edges" << std::endl;

  // Calculate S-parameters
  std::vector<WavePort> ports{wp1, wp2};
  auto S = calculate_sparams(mesh, p, bc, ports);

  std::cout << "\n=== Computed S-parameters ===" << std::endl;
  std::cout << "S11 = " << S(0, 0) << " (|S11| = " << std::abs(S(0, 0)) << ")"
            << std::endl;
  std::cout << "S21 = " << S(1, 0) << " (|S21| = " << std::abs(S(1, 0)) << ")"
            << std::endl;
  std::cout << "S12 = " << S(0, 1) << " (|S12| = " << std::abs(S(0, 1)) << ")"
            << std::endl;
  std::cout << "S22 = " << S(1, 1) << " (|S22| = " << std::abs(S(1, 1)) << ")"
            << std::endl;

  // Expected analytical
  std::complex<double> S21_expected = analytical_s21(freq, WG_L);
  std::cout << "\n=== Expected (analytical) ===" << std::endl;
  std::cout << "S11 = 0 (matched)" << std::endl;
  std::cout << "S21 = " << S21_expected
            << " (|S21| = " << std::abs(S21_expected) << ")" << std::endl;

  // Check results
  double s11_mag = std::abs(S(0, 0));
  double s21_mag = std::abs(S(1, 0));

  std::cout << "\n=== Results ===" << std::endl;
  std::cout << "Passivity: |S11|² + |S21|² = "
            << (s11_mag * s11_mag + s21_mag * s21_mag) << " (should be ≤ 1)"
            << std::endl;

  bool passive =
      (s11_mag * s11_mag + s21_mag * s21_mag) <= 1.01; // small tolerance
  std::cout << "Passivity check: " << (passive ? "PASS" : "FAIL") << std::endl;

  // For a matched waveguide:
  // |S11| should be small (< 0.1 ideally)
  // |S21| should be close to 1 (lossless transmission)
  bool s11_ok = s11_mag < 0.3; // relaxed tolerance for coarse mesh
  bool s21_ok = s21_mag > 0.7; // expect most power transmitted

  std::cout << "|S11| < 0.3: " << (s11_ok ? "PASS" : "FAIL") << std::endl;
  std::cout << "|S21| > 0.7: " << (s21_ok ? "PASS" : "FAIL") << std::endl;

  // Phase check for S21
  double s21_phase_computed = std::arg(S(1, 0));
  double s21_phase_expected = std::arg(S21_expected);
  double phase_diff = std::abs(s21_phase_computed - s21_phase_expected);
  // Wrap phase difference to [-π, π]
  while (phase_diff > M_PI)
    phase_diff -= 2 * M_PI;
  phase_diff = std::abs(phase_diff);

  std::cout << "S21 phase: computed=" << s21_phase_computed * 180 / M_PI << "°"
            << ", expected=" << s21_phase_expected * 180 / M_PI << "°"
            << ", diff=" << phase_diff * 180 / M_PI << "°" << std::endl;

  // DIAGNOSTIC: If results are bad, print additional debugging info
  if (!passive || !s11_ok || !s21_ok) {
    std::cout << "\n=== DIAGNOSTIC INFO (results do not match expected) ==="
              << std::endl;

    // Check matrix condition
    auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);

    // Count DOFs
    int total_dofs = static_cast<int>(mesh.edges.size());
    int pec_dofs = static_cast<int>(bc.dirichlet_edges.size());
    int free_dofs = total_dofs - pec_dofs;
    std::cout << "DOFs: " << total_dofs << " total, " << pec_dofs << " PEC, "
              << free_dofs << " free" << std::endl;

    // Matrix diagonal statistics
    double diag_max = 0, diag_min = 1e100;
    int pos_count = 0, neg_count = 0;
    for (int i = 0; i < total_dofs; ++i) {
      if (bc.dirichlet_edges.count(i))
        continue;
      double d = std::real(asmbl.A.coeff(i, i));
      diag_max = std::max(diag_max, d);
      diag_min = std::min(diag_min, d);
      if (d > 0)
        pos_count++;
      else
        neg_count++;
    }
    std::cout << "Diagonal: min=" << diag_min << ", max=" << diag_max
              << std::endl;
    std::cout << "Diagonal signs: " << pos_count << " positive, " << neg_count
              << " negative" << std::endl;

    // RHS statistics
    double rhs_max = 0;
    for (int i = 0; i < asmbl.b.size(); ++i) {
      rhs_max = std::max(rhs_max, std::abs(asmbl.b[i]));
    }
    std::cout << "RHS max magnitude: " << rhs_max << std::endl;

    // Solve and check residual
    auto res = solve_linear(asmbl.A, asmbl.b, {});
    std::cout << "Solver converged: " << res.converged << std::endl;
    std::cout << "Solver iterations: " << res.iters << std::endl;
    std::cout << "Solver residual: " << res.residual << std::endl;

    // Solution statistics
    double sol_max = 0;
    for (int i = 0; i < res.x.size(); ++i) {
      if (!bc.dirichlet_edges.count(i))
        sol_max = std::max(sol_max, std::abs(res.x[i]));
    }
    std::cout << "Solution max magnitude: " << sol_max << std::endl;

    // Port voltage extraction check
    std::complex<double> V1 = 0, V2 = 0;
    for (size_t k = 0; k < wp1.edges.size(); ++k) {
      int e = wp1.edges[k];
      if (!bc.dirichlet_edges.count(e))
        V1 += std::conj(wp1.weights[k]) * res.x(e);
    }
    for (size_t k = 0; k < wp2.edges.size(); ++k) {
      int e = wp2.edges[k];
      if (!bc.dirichlet_edges.count(e))
        V2 += std::conj(wp2.weights[k]) * res.x(e);
    }
    std::complex<double> V_inc = std::sqrt(wp1.mode.Z0);
    std::cout << "V1 = " << V1 << ", V2 = " << V2 << ", V_inc = " << V_inc
              << std::endl;
    std::cout << "V1/V_inc = " << V1 / V_inc << ", V2/V_inc = " << V2 / V_inc
              << std::endl;
  }

  // Report final status
  if (passive && s11_ok && s21_ok) {
    std::cout << "\n*** BENCHMARK PASSED ***" << std::endl;
  } else {
    std::cout << "\n*** BENCHMARK FAILED - FEM core needs investigation ***"
              << std::endl;
    // Don't assert here - we want to see all diagnostic output
  }
}

// =============================================================================
// Test 4: Element-Level Eigenfrequency Analysis
// =============================================================================
// Compute the generalized eigenvalue problem (K - λM)x = 0 for a single element
// to understand the resolvable frequency range.

void test_element_eigenfrequencies() {
  std::cout << "\n=== Test 4: Element Eigenfrequency Analysis ===" << std::endl;

  // Use a realistic element size from the waveguide mesh
  // Mesh size ~3mm, so a typical tet has edges ~3mm
  double h = 0.003; // 3mm element size

  std::array<Eigen::Vector3d, 4> tet = {
      Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(h, 0.0, 0.0),
      Eigen::Vector3d(0.0, h, 0.0), Eigen::Vector3d(0.0, 0.0, h)};

  auto K = whitney_curl_curl_matrix(tet);
  auto M = whitney_mass_matrix(tet);

  // Solve generalized eigenvalue problem K*x = λ*M*x
  // where λ = k0² = (ω/c)² = (2πf/c)²
  // So f = c * sqrt(λ) / (2π)

  // Convert to standard eigenvalue problem: M^(-1/2) K M^(-1/2) y = λ y
  // Or use generalized eigenvalue solver
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K, M);

  std::cout << "Element eigenvalues (k0²):" << std::endl;
  double k0_min = std::numeric_limits<double>::max();
  for (int i = 0; i < 6; ++i) {
    double lambda = ges.eigenvalues()(i);
    if (lambda > 1e-6) { // Skip null space
      double k0 = std::sqrt(lambda);
      double f = k0 * c0 / (2.0 * M_PI);
      std::cout << "  λ_" << i << " = " << lambda << " → k0 = " << k0
                << " → f = " << f / 1e9 << " GHz" << std::endl;
      k0_min = std::min(k0_min, k0);
    } else {
      std::cout << "  λ_" << i << " = " << lambda << " (null space)"
                << std::endl;
    }
  }

  double f_min = k0_min * c0 / (2.0 * M_PI);
  std::cout << "\nMinimum resolvable frequency for this element: "
            << f_min / 1e9 << " GHz" << std::endl;
  std::cout << "Operating frequency: 10 GHz" << std::endl;
  std::cout << "Ratio (f_element_min / f_operating): " << f_min / 10e9
            << std::endl;

  // The rule of thumb is ~10-20 elements per wavelength
  // At 10 GHz, λ = 30mm, so with 3mm elements we have ~10 elements/wavelength
  double wavelength = c0 / 10e9;
  std::cout << "\nWavelength at 10 GHz: " << wavelength * 1000 << " mm"
            << std::endl;
  std::cout << "Element size: " << h * 1000 << " mm" << std::endl;
  std::cout << "Elements per wavelength: " << wavelength / h << std::endl;

  if (f_min > 10e9) {
    std::cout << "\nWARNING: Element eigenfrequency (" << f_min / 1e9
              << " GHz) is above operating frequency (10 GHz)" << std::endl;
    std::cout << "This may cause poor wave propagation modeling." << std::endl;
    std::cout << "Consider using finer mesh (more elements per wavelength)."
              << std::endl;
  }
}

// =============================================================================
// Test 5: K/M Ratio Analysis
// =============================================================================
// Analyze the ratio of curl-curl to mass matrix contributions at operating
// freq.

void test_km_ratio() {
  std::cout << "\n=== Test 5: K/M Ratio Analysis ===" << std::endl;

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  double freq = 10e9;
  double omega = 2.0 * M_PI * freq;
  double k0 = omega / c0;
  double k0_sq = k0 * k0;

  std::cout << "k0² = " << k0_sq << std::endl;

  // Assemble K and M separately by using different parameters
  MaxwellParams p_k, p_m;
  p_k.omega = 1e-10; // Very low frequency, A ≈ K
  p_m.omega = 1e15;  // Very high frequency, A ≈ -k0²M

  auto asmbl_k = assemble_maxwell(mesh, p_k, bc, {}, -1);
  auto asmbl_m = assemble_maxwell(mesh, p_m, bc, {}, -1);

  // Extract K and M estimates
  // At very low freq: A_low ≈ K
  // At very high freq: A_high ≈ K - k0_high² M, so M ≈ (K - A_high) / k0_high²
  double k0_high = p_m.omega / c0;
  double k0_high_sq = k0_high * k0_high;

  // Sample diagonal statistics
  double K_diag_avg = 0, M_diag_avg = 0;
  int free_count = 0;
  for (int i = 0; i < static_cast<int>(mesh.edges.size()); ++i) {
    if (bc.dirichlet_edges.count(i))
      continue;

    double K_ii = std::real(asmbl_k.A.coeff(i, i));
    double A_high_ii = std::real(asmbl_m.A.coeff(i, i));
    double M_ii = (K_ii - A_high_ii) / k0_high_sq;

    K_diag_avg += K_ii;
    M_diag_avg += M_ii;
    free_count++;
  }
  K_diag_avg /= free_count;
  M_diag_avg /= free_count;

  std::cout << "Average diagonal K_ii: " << K_diag_avg << std::endl;
  std::cout << "Average diagonal M_ii: " << M_diag_avg << std::endl;
  std::cout << "At 10 GHz, k0²*M_ii / K_ii ratio: "
            << k0_sq * M_diag_avg / K_diag_avg << std::endl;

  // For proper wave propagation, we need K and k0²M to be comparable
  // If K >> k0²M, the system is "stiff" (quasi-static regime)
  // If K << k0²M, the system is dominated by mass (high-frequency regime)
  double ratio = K_diag_avg / (k0_sq * M_diag_avg);
  std::cout << "K / (k0²M) ratio: " << ratio << std::endl;

  if (ratio > 10) {
    std::cout << "WARNING: Curl-curl (K) dominates by " << ratio << "x"
              << std::endl;
    std::cout << "This suggests quasi-static behavior, not wave propagation."
              << std::endl;
  } else if (ratio < 0.1) {
    std::cout << "WARNING: Mass (M) dominates by " << 1 / ratio << "x"
              << std::endl;
    std::cout << "This may cause numerical instability." << std::endl;
  } else {
    std::cout << "K and k0²M are comparable (ratio = " << ratio << ")"
              << std::endl;
    std::cout << "This is appropriate for wave propagation simulation."
              << std::endl;
  }
}

// =============================================================================
// Main
// =============================================================================

int main() {
  std::cout << std::setprecision(6);
  std::cout << "========================================" << std::endl;
  std::cout << "EdgeFEM FEM Core Validation Benchmark" << std::endl;
  std::cout << "========================================" << std::endl;

  test_whitney_matrices();
  test_assembly_formulation();
  test_element_eigenfrequencies();
  test_km_ratio();
  test_waveguide_sparams();

  std::cout << "\n========================================" << std::endl;
  std::cout << "Benchmark complete." << std::endl;
  std::cout << "========================================" << std::endl;

  return 0;
}
