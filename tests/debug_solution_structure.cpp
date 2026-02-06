// Examine the FEM solution structure when exciting with port source
// Key question: Does the solution look like a traveling wave?

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <Eigen/Dense>
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
  double beta = std::sqrt(k0*k0 - std::pow(M_PI/0.02286, 2));  // Propagation constant

  RectWaveguidePort dims{0.02286, 0.01016};

  PortSurfaceMesh port1_surf = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2_surf = extract_surface_mesh(mesh, 3);

  PortMode mode1 = solve_te10_mode(dims, freq);
  PortMode mode2 = solve_te10_mode(dims, freq);
  populate_te10_field(port1_surf, dims, mode1);
  populate_te10_field(port2_surf, dims, mode2);

  WavePort wp1 = build_wave_port(mesh, port1_surf, mode1);
  WavePort wp2 = build_wave_port(mesh, port2_surf, mode2);

  std::cout << "=== Waveguide Parameters ===" << std::endl;
  std::cout << "Frequency: " << freq/1e9 << " GHz" << std::endl;
  std::cout << "k0 = " << k0 << " rad/m" << std::endl;
  std::cout << "beta = " << beta << " rad/m (expected: " << std::real(mode1.beta) << ")" << std::endl;
  std::cout << "Wavelength = " << 2*M_PI/beta*1000 << " mm" << std::endl;
  std::cout << "Waveguide length L = 50 mm" << std::endl;
  std::cout << "Phase shift beta*L = " << beta*0.05 << " rad = " << beta*0.05*180/M_PI << " deg" << std::endl;

  // Solve with ports
  MaxwellParams p;
  p.omega = omega;

  std::vector<WavePort> ports{wp1, wp2};
  auto asmbl = assemble_maxwell(mesh, p, bc, ports, 0);  // Excite port 1
  auto res = solve_linear(asmbl.A, asmbl.b, {});

  std::cout << "\n=== Solution Statistics ===" << std::endl;
  std::cout << "||x|| = " << res.x.norm() << std::endl;
  std::cout << "||b|| = " << asmbl.b.norm() << std::endl;

  // Analyze solution magnitude vs z
  std::map<double, std::vector<std::pair<int, std::complex<double>>>> z_groups;

  for (int e = 0; e < (int)mesh.edges.size(); ++e) {
    if (bc.dirichlet_edges.count(e)) continue;

    const auto &edge = mesh.edges[e];
    const auto &p0 = mesh.nodes.at(mesh.nodeIndex.at(edge.n0)).xyz;
    const auto &p1 = mesh.nodes.at(mesh.nodeIndex.at(edge.n1)).xyz;
    double z_mid = (p0.z() + p1.z()) / 2.0;
    double z_rounded = std::round(z_mid * 1000) / 1000.0;

    z_groups[z_rounded].push_back({e, res.x(e)});
  }

  std::cout << "\n=== Solution Magnitude vs Z ===" << std::endl;
  std::cout << "Z(mm)\t\tNum Edges\tMean |x|\t\tPhase (deg)" << std::endl;

  // For traveling wave exp(-j*beta*z), phase should increase linearly with z
  double prev_phase = 0;
  bool first = true;

  for (const auto &kv : z_groups) {
    double z = kv.first;
    const auto &edges = kv.second;

    // Compute mean magnitude and phase
    double sum_mag = 0;
    std::complex<double> sum_x = 0;
    for (const auto &ex : edges) {
      sum_mag += std::abs(ex.second);
      sum_x += ex.second;
    }
    double mean_mag = sum_mag / edges.size();
    double phase = std::arg(sum_x) * 180 / M_PI;

    // Only show some z-values
    if (z < 0.003 || z > 0.047 || (z > 0.022 && z < 0.028)) {
      std::cout << z*1000 << "\t\t" << edges.size() << "\t\t" << mean_mag << "\t\t" << phase << std::endl;
    }

    if (first) {
      prev_phase = phase;
      first = false;
    }
  }

  // Compare port 1 and port 2 solutions in detail
  std::cout << "\n=== Port 1 Solution Detail ===" << std::endl;
  std::complex<double> sum_x1 = 0, sum_w1x1 = 0;
  int count1 = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;
    sum_x1 += res.x(e);
    sum_w1x1 += std::conj(wp1.weights(i)) * res.x(e);
    count1++;
  }
  std::cout << "Mean x at port 1: " << sum_x1 / (double)count1 << std::endl;
  std::cout << "V1 = w^H * x = " << sum_w1x1 << std::endl;

  std::cout << "\n=== Port 2 Solution Detail ===" << std::endl;
  std::complex<double> sum_x2 = 0, sum_w2x2 = 0;
  int count2 = 0;
  for (size_t i = 0; i < wp2.edges.size(); ++i) {
    int e = wp2.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;
    sum_x2 += res.x(e);
    sum_w2x2 += std::conj(wp2.weights(i)) * res.x(e);
    count2++;
  }
  std::cout << "Mean x at port 2: " << sum_x2 / (double)count2 << std::endl;
  std::cout << "V2 = w^H * x = " << sum_w2x2 << std::endl;

  // Check if solution at port 2 relates to port 1 by exp(-j*beta*L)
  std::complex<double> expected_ratio = std::exp(std::complex<double>(0, -beta * 0.05));
  std::cout << "\nExpected ratio x2/x1 for traveling wave: " << expected_ratio << std::endl;
  std::cout << "Actual ratio (mean x): " << sum_x2/sum_x1 << std::endl;

  // Analyze the source term
  std::cout << "\n=== Source Term Analysis ===" << std::endl;
  int nonzero_b = 0;
  double b_sum = 0;
  for (int e = 0; e < (int)mesh.edges.size(); ++e) {
    if (std::abs(asmbl.b(e)) > 1e-15) {
      nonzero_b++;
      b_sum += std::abs(asmbl.b(e));
    }
  }
  std::cout << "Non-zero source entries: " << nonzero_b << std::endl;
  std::cout << "Port 1 free edges: " << (wp1.edges.size() - 24) << std::endl;
  std::cout << "Total source magnitude: " << b_sum << std::endl;

  // Check the matrix structure at port edges
  std::cout << "\n=== Matrix Structure at Ports ===" << std::endl;

  // Find a sample port edge
  int sample_edge = -1;
  for (int e : wp1.edges) {
    if (!bc.dirichlet_edges.count(e)) {
      sample_edge = e;
      break;
    }
  }

  if (sample_edge >= 0) {
    std::cout << "Sample port edge: " << sample_edge << std::endl;

    // Count non-zero entries in row
    int nnz_row = 0;
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A, sample_edge); it; ++it) {
      if (std::abs(it.value()) > 1e-15) nnz_row++;
    }
    std::cout << "Non-zeros in row: " << nnz_row << std::endl;

    // Diagonal value
    std::cout << "Diagonal A(" << sample_edge << "," << sample_edge << ") = " << asmbl.A.coeff(sample_edge, sample_edge) << std::endl;

    // Source value
    std::cout << "Source b(" << sample_edge << ") = " << asmbl.b(sample_edge) << std::endl;

    // Solution value
    std::cout << "Solution x(" << sample_edge << ") = " << res.x(sample_edge) << std::endl;
  }

  // Key insight: compare the port loading contribution
  std::cout << "\n=== Port Loading Contribution ===" << std::endl;
  double Z0 = std::real(wp1.mode.Z0);
  std::cout << "Z0 = " << Z0 << " ohms" << std::endl;
  std::cout << "||w||^2 = " << wp1.weights.squaredNorm() << std::endl;
  std::cout << "||w||^2 / Z0 = " << wp1.weights.squaredNorm() / Z0 << std::endl;

  // The port loading adds ww^H/Z0 to the matrix
  // For edge e in port: A(e,e) += |w_e|^2 / Z0
  // This should be comparable to the original A(e,e) for proper coupling

  double sum_w_sq = 0, sum_diag = 0;
  for (size_t i = 0; i < wp1.edges.size(); ++i) {
    int e = wp1.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;
    double w_sq = std::norm(wp1.weights(i));
    sum_w_sq += w_sq;

    // Original diagonal (A without port loading)
    // We can't easily get this, so skip
  }
  std::cout << "Sum |w_e|^2 for port edges = " << sum_w_sq << std::endl;
  std::cout << "Ratio to ||w||^2: " << sum_w_sq / wp1.weights.squaredNorm() << std::endl;

  return 0;
}
