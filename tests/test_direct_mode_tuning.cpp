// Tune the direct eigenmode excitation for better S-parameter accuracy
// Try different ABC coefficients and normalization approaches

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_set>

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);
constexpr double eta0 = mu0 * c0;

struct TestResult {
  std::string name;
  double abc_scale;
  std::complex<double> S11, S21;
  double passivity;
};

void run_test(const Mesh &mesh, const BC &bc, const PortSurfaceMesh &port1_surf,
              const PortSurfaceMesh &port2_surf, const Eigen::VectorXcd &v1,
              const Eigen::VectorXcd &v2, double omega, double beta,
              double abc_scale, const std::string &name,
              std::vector<TestResult> &results) {

  int m = mesh.edges.size();

  // Collect port edges
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

  // Assemble Maxwell system
  MaxwellParams p;
  p.omega = omega;

  std::vector<WavePort> no_ports;
  auto asmbl = assemble_maxwell(mesh, p, bc, no_ports, -1);

  // Add ABC with scaled coefficient
  const std::complex<double> abc_coeff(0.0, beta * abc_scale);

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

  // Excitation: RHS = 2*j*beta*abc_scale * v1
  asmbl.b = 2.0 * abc_coeff * v1;

  // Solve
  auto res = solve_linear(asmbl.A, asmbl.b, {});

  // S-parameters
  std::complex<double> V1 = v1.dot(res.x);
  std::complex<double> V2 = v2.dot(res.x);
  std::complex<double> V_inc(1.0, 0.0);

  std::complex<double> S11 = (V1 - V_inc) / V_inc;
  std::complex<double> S21 = V2 / V_inc;

  TestResult r;
  r.name = name;
  r.abc_scale = abc_scale;
  r.S11 = S11;
  r.S21 = S21;
  r.passivity = std::norm(S11) + std::norm(S21);
  results.push_back(r);
}

int main() {
  std::cout << std::setprecision(5) << std::fixed;

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);
  int m = mesh.edges.size();

  double freq = 10e9;
  double omega = 2 * M_PI * freq;
  double k0 = omega / c0;
  double L = 0.05;

  RectWaveguidePort dims{0.02286, 0.01016};
  double kc = M_PI / dims.a;
  double beta = std::sqrt(k0 * k0 - kc * kc);
  double Z0_te10 = omega * mu0 / beta;

  std::cout << "=== ABC Coefficient Tuning Test ===" << std::endl;
  std::cout << "Target: |S11| < 0.1, |S21| > 0.95" << std::endl;
  std::cout << "Expected phase(S21) = " << -beta * L * 180 / M_PI << "°"
            << std::endl;
  std::cout << std::endl;

  // Extract port surfaces
  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  // Compute eigenvector
  double kc_sq = kc * kc;
  Eigen::VectorXd v_te10 =
      compute_te_eigenvector(mesh, bc.dirichlet_edges, kc_sq);

  // Create port eigenvectors (normalized)
  Eigen::VectorXcd v1 = Eigen::VectorXcd::Zero(m);
  Eigen::VectorXcd v2 = Eigen::VectorXcd::Zero(m);

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

  for (int e : port1_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      v1(e) = std::complex<double>(0.0, v_te10(e));
    }
  }
  for (int e : port2_edges) {
    if (!bc.dirichlet_edges.count(e)) {
      v2(e) = std::complex<double>(0.0, v_te10(e));
    }
  }

  double v1_norm = v1.norm();
  double v2_norm = v2.norm();
  if (v1_norm > 1e-10)
    v1 /= v1_norm;
  if (v2_norm > 1e-10)
    v2 /= v2_norm;

  std::vector<TestResult> results;

  // Try different ABC scaling factors
  for (double scale :
       {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0}) {
    std::string name = "ABC scale=" + std::to_string(scale);
    run_test(mesh, bc, port1_surf, port2_surf, v1, v2, omega, beta, scale, name,
             results);
  }

  // Try with frequency-normalized ABC: j*beta/k0
  for (double scale : {0.5, 1.0, 2.0}) {
    double abc_eff = beta / k0 * scale;
    std::string name = "ABC beta/k0 scale=" + std::to_string(scale);
    run_test(mesh, bc, port1_surf, port2_surf, v1, v2, omega, beta,
             abc_eff / beta, name, results);
  }

  // Try impedance-matched ABC: j*beta*sqrt(Z0/eta0)
  double imp_scale = std::sqrt(Z0_te10 / eta0);
  for (double scale : {0.5, 1.0, 2.0}) {
    std::string name = "ABC imp-match scale=" + std::to_string(scale);
    run_test(mesh, bc, port1_surf, port2_surf, v1, v2, omega, beta,
             imp_scale * scale, name, results);
  }

  // Print results
  std::cout << std::setw(30) << std::left << "Configuration" << std::setw(12)
            << "|S11|" << std::setw(12) << "|S21|" << std::setw(12)
            << "phase(S21)" << std::setw(12) << "passivity" << std::endl;
  std::cout << std::string(78, '-') << std::endl;

  for (const auto &r : results) {
    double phase_deg = std::arg(r.S21) * 180.0 / M_PI;
    std::cout << std::setw(30) << std::left << r.name << std::setw(12)
              << std::abs(r.S11) << std::setw(12) << std::abs(r.S21)
              << std::setw(12) << phase_deg << std::setw(12) << r.passivity;
    if (std::abs(r.S11) < 0.1 && std::abs(r.S21) > 0.95 &&
        r.passivity <= 1.01) {
      std::cout << " ✓ PASS";
    }
    std::cout << std::endl;
  }

  // Find best configuration
  double best_metric = 0;
  std::string best_name;
  for (const auto &r : results) {
    // Metric: maximize |S21| while minimizing |S11|, penalize passivity
    // violation
    double metric = std::abs(r.S21) - std::abs(r.S11);
    if (r.passivity > 1.01)
      metric -= (r.passivity - 1.0);
    if (metric > best_metric) {
      best_metric = metric;
      best_name = r.name;
    }
  }
  std::cout << "\nBest configuration: " << best_name << std::endl;

  return 0;
}
