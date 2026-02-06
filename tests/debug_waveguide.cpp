#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <set>

#include "edgefem/maxwell.hpp"
#include "edgefem/solver.hpp"

using namespace edgefem;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  std::cout << "Mesh loaded: " << mesh.nodes.size() << " nodes, "
            << mesh.edges.size() << " edges, "
            << mesh.tets.size() << " tets" << std::endl;

  BC bc = build_edge_pec(mesh, 1);
  std::cout << "PEC BC: " << bc.dirichlet_edges.size() << " edges" << std::endl;

  const double freq = 10e9;
  MaxwellParams p;
  p.omega = 2 * M_PI * freq;
  std::cout << "Frequency: " << freq/1e9 << " GHz, omega = " << p.omega << std::endl;

  PortSurfaceMesh port1 = extract_surface_mesh(mesh, 2);
  PortSurfaceMesh port2 = extract_surface_mesh(mesh, 3);
  std::cout << "Port 1: " << port1.mesh.nodes.size() << " nodes, "
            << port1.mesh.tris.size() << " tris" << std::endl;
  std::cout << "Port 2: " << port2.mesh.nodes.size() << " nodes, "
            << port2.mesh.tris.size() << " tris" << std::endl;

  // Eigensolve-based modes
  std::cout << "\n=== Eigensolve modes ===" << std::endl;
  auto modes1 = solve_port_eigens(port1.mesh, 1, p.omega, 1.0, 1.0, ModePolarization::TE);
  auto modes2 = solve_port_eigens(port2.mesh, 1, p.omega, 1.0, 1.0, ModePolarization::TE);

  if (!modes1.empty()) {
    const auto& m = modes1.front();
    std::cout << "Mode 1: fc=" << m.fc/1e9 << " GHz, kc=" << m.kc
              << ", beta=" << m.beta << ", Z0=" << m.Z0 << std::endl;
  }
  if (!modes2.empty()) {
    const auto& m = modes2.front();
    std::cout << "Mode 2: fc=" << m.fc/1e9 << " GHz, kc=" << m.kc
              << ", beta=" << m.beta << ", Z0=" << m.Z0 << std::endl;
  }

  // Analytical modes
  std::cout << "\n=== Analytical modes ===" << std::endl;
  RectWaveguidePort dims{0.02286, 0.01016};
  PortMode mode1_ana = solve_te10_mode(dims, freq);
  PortMode mode2_ana = solve_te10_mode(dims, freq);
  std::cout << "Analytical: fc=" << mode1_ana.fc/1e9 << " GHz, kc=" << mode1_ana.kc
            << ", beta=" << mode1_ana.beta << ", Z0=" << mode1_ana.Z0 << std::endl;

  // Build wave ports using eigensolve modes
  if (!modes1.empty() && !modes2.empty()) {
    WavePort wp1 = build_wave_port(mesh, port1, modes1.front());
    WavePort wp2 = build_wave_port(mesh, port2, modes2.front());
    std::cout << "\nWave port 1: " << wp1.edges.size() << " edges" << std::endl;
    std::cout << "Wave port 2: " << wp2.edges.size() << " edges" << std::endl;

    double max_w1 = 0, max_w2 = 0;
    for (int i = 0; i < wp1.weights.size(); ++i)
      max_w1 = std::max(max_w1, std::abs(wp1.weights[i]));
    for (int i = 0; i < wp2.weights.size(); ++i)
      max_w2 = std::max(max_w2, std::abs(wp2.weights[i]));
    std::cout << "Max |weight| port1: " << max_w1 << std::endl;
    std::cout << "Max |weight| port2: " << max_w2 << std::endl;

    std::vector<WavePort> ports{wp1, wp2};
    auto S = calculate_sparams(mesh, p, bc, ports);

    std::cout << "\n=== S-parameters (eigensolve modes) ===" << std::endl;
    std::cout << "S11 = " << S(0,0) << " (|S11| = " << std::abs(S(0,0)) << ")" << std::endl;
    std::cout << "S21 = " << S(1,0) << " (|S21| = " << std::abs(S(1,0)) << ")" << std::endl;
    std::cout << "S12 = " << S(0,1) << " (|S12| = " << std::abs(S(0,1)) << ")" << std::endl;
    std::cout << "S22 = " << S(1,1) << " (|S22| = " << std::abs(S(1,1)) << ")" << std::endl;
  }

  // Build wave ports using analytical modes
  std::cout << "\n=== Analytical mode S-parameters ===" << std::endl;
  populate_te10_field(port1, dims, mode1_ana);
  populate_te10_field(port2, dims, mode2_ana);

  // Verify mode field values
  std::cout << "\n=== Mode field verification ===" << std::endl;
  double x_min = 1e10, x_max = -1e10;
  for (const auto& node : port1.mesh.nodes) {
    x_min = std::min(x_min, node.xyz.x());
    x_max = std::max(x_max, node.xyz.x());
  }
  std::cout << "Port 1 x-range: [" << x_min*1000 << ", " << x_max*1000 << "] mm" << std::endl;
  std::cout << "Sample Hz values vs analytical:" << std::endl;
  for (int i = 0; i < std::min(5, (int)port1.mesh.nodes.size()); ++i) {
    const auto& node = port1.mesh.nodes[i];
    double x = node.xyz.x() - x_min;
    double Hz_expected = std::cos(M_PI * x / dims.a);
    std::cout << "  Node " << i << ": x=" << node.xyz.x()*1000 << "mm"
              << ", Hz=" << mode1_ana.field(i)
              << ", expected=" << Hz_expected << std::endl;
  }

  // Check port surface z-coordinates
  double port1_z_min = 1e10, port1_z_max = -1e10;
  double port2_z_min = 1e10, port2_z_max = -1e10;
  for (int idx : port1.volume_tri_indices) {
    const auto& tri = mesh.tris[idx];
    for (int k = 0; k < 3; ++k) {
      double z = mesh.nodes[mesh.nodeIndex.at(tri.conn[k])].xyz.z();
      port1_z_min = std::min(port1_z_min, z);
      port1_z_max = std::max(port1_z_max, z);
    }
  }
  for (int idx : port2.volume_tri_indices) {
    const auto& tri = mesh.tris[idx];
    for (int k = 0; k < 3; ++k) {
      double z = mesh.nodes[mesh.nodeIndex.at(tri.conn[k])].xyz.z();
      port2_z_min = std::min(port2_z_min, z);
      port2_z_max = std::max(port2_z_max, z);
    }
  }
  std::cout << "Port 1 z-range: [" << port1_z_min << ", " << port1_z_max << "]" << std::endl;
  std::cout << "Port 2 z-range: [" << port2_z_min << ", " << port2_z_max << "]" << std::endl;

  WavePort wp1_ana = build_wave_port(mesh, port1, mode1_ana);
  WavePort wp2_ana = build_wave_port(mesh, port2, mode2_ana);

  // Test different weight scaling factors
  std::cout << "\n=== Weight scaling experiment ===" << std::endl;
  std::vector<double> scale_factors = {1.0, 0.1, 0.01, 0.001};
  for (double sf : scale_factors) {
    WavePort wp1_test = wp1_ana;
    WavePort wp2_test = wp2_ana;
    wp1_test.weights *= sf;
    wp2_test.weights *= sf;

    std::vector<WavePort> ports_test{wp1_test, wp2_test};
    auto S_test = calculate_sparams(mesh, p, bc, ports_test);

    std::cout << "Scale " << sf << ": S11=" << S_test(0,0)
              << ", S21=" << S_test(1,0) << std::endl;
  }
  std::cout << "=== Manual weight scaling test ===" << std::endl;

  double max_w1_ana = 0, max_w2_ana = 0;
  std::complex<double> sum_w1 = 0, sum_w2 = 0;
  for (int i = 0; i < wp1_ana.weights.size(); ++i) {
    max_w1_ana = std::max(max_w1_ana, std::abs(wp1_ana.weights[i]));
    sum_w1 += wp1_ana.weights[i];
  }
  for (int i = 0; i < wp2_ana.weights.size(); ++i) {
    max_w2_ana = std::max(max_w2_ana, std::abs(wp2_ana.weights[i]));
    sum_w2 += wp2_ana.weights[i];
  }
  std::cout << "Max |weight| port1 (ana): " << max_w1_ana << std::endl;
  std::cout << "Max |weight| port2 (ana): " << max_w2_ana << std::endl;
  std::cout << "Sum of weights port1: " << sum_w1 << std::endl;
  std::cout << "Sum of weights port2: " << sum_w2 << std::endl;
  // Show a few sample weights
  std::cout << "Sample weights port1: ";
  for (int i = 0; i < std::min(5, (int)wp1_ana.weights.size()); ++i)
    std::cout << wp1_ana.weights[i] << " ";
  std::cout << std::endl;

  std::vector<WavePort> ports_ana{wp1_ana, wp2_ana};
  auto S_ana = calculate_sparams(mesh, p, bc, ports_ana);

  std::cout << "S11 = " << S_ana(0,0) << " (|S11| = " << std::abs(S_ana(0,0)) << ")" << std::endl;
  std::cout << "S21 = " << S_ana(1,0) << " (|S21| = " << std::abs(S_ana(1,0)) << ")" << std::endl;
  std::cout << "S12 = " << S_ana(0,1) << " (|S12| = " << std::abs(S_ana(0,1)) << ")" << std::endl;
  std::cout << "S22 = " << S_ana(1,1) << " (|S22| = " << std::abs(S_ana(1,1)) << ")" << std::endl;

  // Expected analytical
  double length = 0.05;
  auto analytic = straight_waveguide_sparams(dims, length, freq);
  std::cout << "\n=== Expected (analytical) ===" << std::endl;
  std::cout << "S11 = " << analytic.s11 << std::endl;
  std::cout << "S21 = " << analytic.s21 << " (|S21| = " << std::abs(analytic.s21) << ")" << std::endl;

  // Debug: check assembly diagnostics
  std::cout << "\n=== Assembly diagnostics ===" << std::endl;
  auto asmbl = assemble_maxwell(mesh, p, bc, ports_ana, 0);

  // Check matrix properties
  double max_diag = 0, min_diag = 1e100;
  std::complex<double> diag_sum = 0;
  int free_dof_count = 0;
  for (int i = 0; i < asmbl.A.rows(); ++i) {
    std::complex<double> d = asmbl.A.coeff(i, i);
    double mag = std::abs(d);
    if (!bc.dirichlet_edges.count(i)) {
      max_diag = std::max(max_diag, mag);
      min_diag = std::min(min_diag, mag);
      diag_sum += d;
      free_dof_count++;
    }
  }
  std::cout << "Free DOFs: " << free_dof_count << std::endl;
  std::cout << "Matrix diagonal range: [" << min_diag << ", " << max_diag << "]" << std::endl;

  // Check RHS
  double max_b = 0;
  int nonzero_b = 0;
  std::complex<double> b_sum = 0;
  for (int i = 0; i < asmbl.b.size(); ++i) {
    double b_mag = std::abs(asmbl.b[i]);
    max_b = std::max(max_b, b_mag);
    b_sum += asmbl.b[i];
    if (b_mag > 1e-15) nonzero_b++;
  }
  std::cout << "Max |b|: " << max_b << ", nonzero entries: " << nonzero_b << std::endl;
  std::cout << "Sum of b: " << b_sum << std::endl;

  // Expected: b = 2*w/sqrt(Z0), so sum(b) = 2*sum(w)/sqrt(Z0)
  std::cout << "Expected sum(b) = 2*sum(w)/sqrt(Z0) = "
            << 2.0 * sum_w1 / std::sqrt(wp1_ana.mode.Z0) << std::endl;

  // Solve and check solution
  auto res = solve_linear(asmbl.A, asmbl.b, {});
  std::cout << "Solver converged: " << res.converged << std::endl;

  // Check solution statistics
  double max_x = 0, sum_x = 0;
  int port1_edge_count = 0;
  std::complex<double> V1_check = 0;
  for (size_t k = 0; k < wp1_ana.edges.size(); ++k) {
    int edge_idx = wp1_ana.edges[k];
    std::complex<double> x_val = res.x(edge_idx);
    std::complex<double> w_val = wp1_ana.weights[k];
    V1_check += std::conj(w_val) * x_val;
    max_x = std::max(max_x, std::abs(x_val));
    sum_x += std::abs(x_val);
    port1_edge_count++;
  }
  std::cout << "\nPort 1 solution check:" << std::endl;
  std::cout << "  Edges on port: " << port1_edge_count << std::endl;
  std::cout << "  Max |x| at port: " << max_x << std::endl;
  std::cout << "  Avg |x| at port: " << sum_x/port1_edge_count << std::endl;
  std::cout << "  V1 = sum(conj(w)*x) = " << V1_check << std::endl;
  std::cout << "  V_inc = sqrt(Z0) = " << std::sqrt(wp1_ana.mode.Z0) << std::endl;
  std::cout << "  V1 / V_inc = " << V1_check / std::sqrt(wp1_ana.mode.Z0) << std::endl;

  // Check global solution
  double global_max_x = 0;
  double global_sum_x = 0;
  for (int i = 0; i < res.x.size(); ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      global_max_x = std::max(global_max_x, std::abs(res.x(i)));
      global_sum_x += std::abs(res.x(i));
    }
  }
  std::cout << "\nGlobal solution:" << std::endl;
  std::cout << "  Max |x|: " << global_max_x << std::endl;
  std::cout << "  Avg |x| (free DOFs): " << global_sum_x/free_dof_count << std::endl;

  // Check port admittance contribution
  std::cout << "\nPort admittance check:" << std::endl;
  double Y_diag_sum = 0;
  for (size_t a = 0; a < wp1_ana.edges.size(); ++a) {
    Y_diag_sum += std::abs(wp1_ana.weights[a] * std::conj(wp1_ana.weights[a]) / wp1_ana.mode.Z0);
  }
  std::cout << "  Sum of Y diagonal: " << Y_diag_sum << std::endl;
  std::cout << "  |w|^2 / Z0 (max weight): " << 390.989*390.989/498.974 << std::endl;

  // Check ratio of port admittance to matrix diagonal
  double avg_diag = (min_diag + max_diag) / 2;
  std::cout << "  Ratio Y_sum / avg_diag: " << Y_diag_sum / avg_diag << std::endl;

  // The relationship for matched port (weak port loading approximation):
  // V = w^H A^-1 b ≈ 2 w^H A^-1 w / sqrt(Z0)
  // For V = V_inc = sqrt(Z0), we need w^H A^-1 w = Z0/2
  // Let's estimate w^H A^-1 w using diagonal approximation
  double wAinvw_est_all = 0, wAinvw_est_free = 0;
  double min_port_diag = 1e20, max_port_diag = 0;
  double min_free_diag = 1e20, max_free_diag = 0;
  double sum_wsq = 0, sum_wsq_free = 0;
  int pec_count = 0, free_count = 0;
  for (size_t a = 0; a < wp1_ana.edges.size(); ++a) {
    int ga = wp1_ana.edges[a];
    std::complex<double> wa = wp1_ana.weights[a];
    std::complex<double> Aga = asmbl.A.coeff(ga, ga);
    double wsq = std::real(std::conj(wa) * wa);
    sum_wsq += wsq;
    double Aga_mag = std::abs(Aga);
    min_port_diag = std::min(min_port_diag, Aga_mag);
    max_port_diag = std::max(max_port_diag, Aga_mag);
    if (Aga_mag > 1e-10) {
      wAinvw_est_all += wsq / Aga_mag;
    }
    // Check if PEC edge
    if (bc.dirichlet_edges.count(ga)) {
      pec_count++;
    } else {
      free_count++;
      sum_wsq_free += wsq;
      min_free_diag = std::min(min_free_diag, Aga_mag);
      max_free_diag = std::max(max_free_diag, Aga_mag);
      if (Aga_mag > 1e-10) {
        wAinvw_est_free += wsq / Aga_mag;
      }
    }
  }
  std::cout << "  Port edges: " << wp1_ana.edges.size() << " total, "
            << pec_count << " PEC, " << free_count << " free" << std::endl;
  std::cout << "  All port edge diagonal range: [" << min_port_diag << ", " << max_port_diag << "]" << std::endl;
  std::cout << "  Free port edge diagonal range: [" << min_free_diag << ", " << max_free_diag << "]" << std::endl;
  std::cout << "  Sum |w|^2 (all): " << sum_wsq << std::endl;
  std::cout << "  Sum |w|^2 (free): " << sum_wsq_free << std::endl;
  std::cout << "  Estimated w^H A^-1 w (all, diag approx): " << wAinvw_est_all << std::endl;
  std::cout << "  Estimated w^H A^-1 w (free, diag approx): " << wAinvw_est_free << std::endl;
  std::cout << "  Target for matching (Z0/2): " << std::real(wp1_ana.mode.Z0)/2 << std::endl;
  std::cout << "  Ratio (free/target): " << wAinvw_est_free / (std::real(wp1_ana.mode.Z0)/2) << std::endl;

  // Check a few PEC edges
  int pec_sample = 0;
  for (int e : bc.dirichlet_edges) {
    if (pec_sample++ < 3) {
      std::cout << "PEC edge " << e << ": A[e,e]=" << asmbl.A.coeff(e,e)
                << ", b[e]=" << asmbl.b[e] << ", x[e]=" << res.x(e) << std::endl;
    }
  }

  // Verify port loading is applied by checking matrix without ports
  std::cout << "\n=== Port loading verification ===" << std::endl;
  auto asmbl_no_port = assemble_maxwell(mesh, p, bc, {}, -1);

  // Check if diagonal changed at port edges
  double max_diag_change = 0;
  for (int ga : wp1_ana.edges) {
    if (!bc.dirichlet_edges.count(ga)) {
      double A_with = std::abs(asmbl.A.coeff(ga, ga));
      double A_without = std::abs(asmbl_no_port.A.coeff(ga, ga));
      max_diag_change = std::max(max_diag_change, std::abs(A_with - A_without));
    }
  }
  std::cout << "Max diagonal change at port edges (Y contribution): " << max_diag_change << std::endl;

  // Check expected Y diagonal contribution
  double expected_Y_diag = 0;
  for (int i = 0; i < wp1_ana.weights.size(); ++i) {
    if (!bc.dirichlet_edges.count(wp1_ana.edges[i])) {
      std::complex<double> w = wp1_ana.weights[i];
      expected_Y_diag = std::max(expected_Y_diag,
                                 std::abs(w * std::conj(w) / wp1_ana.mode.Z0));
    }
  }
  std::cout << "Max expected Y diagonal: " << expected_Y_diag << std::endl;

  // Check solution at port 2
  std::cout << "\n=== Port 2 solution check ===" << std::endl;
  std::complex<double> V2_check = 0;
  for (size_t k = 0; k < wp2_ana.edges.size(); ++k) {
    int edge_idx = wp2_ana.edges[k];
    if (!bc.dirichlet_edges.count(edge_idx)) {
      V2_check += std::conj(wp2_ana.weights[k]) * res.x(edge_idx);
    }
  }
  std::cout << "V2 = sum(conj(w)*x) = " << V2_check << std::endl;
  std::cout << "Expected S21 from V2: " << V2_check / std::sqrt(wp1_ana.mode.Z0) << std::endl;

  // Check edge solution values at both ports
  std::cout << "\n=== Sample edge solutions ===" << std::endl;
  std::cout << "Port 1 (first 10 free edges):" << std::endl;
  int shown = 0;
  for (int i = 0; i < (int)wp1_ana.edges.size() && shown < 10; ++i) {
    int e = wp1_ana.edges[i];
    if (bc.dirichlet_edges.count(e)) continue;
    const auto& edge = mesh.edges[e];
    const auto& pa = mesh.nodes[mesh.nodeIndex.at(edge.n0)].xyz;
    const auto& pb = mesh.nodes[mesh.nodeIndex.at(edge.n1)].xyz;
    Eigen::Vector3d midpt = (pa + pb) / 2.0;
    // For TE10: Ey = sin(pi*x/a), so weight should be ~sin(pi*x/a)*edge_y_length
    double x_norm = midpt.x() / 0.02286; // x/a
    double expected_pattern = std::sin(M_PI * x_norm);
    double edge_y_len = std::abs(pb.y() - pa.y());
    std::cout << "  Edge " << e << ": midpt=(" << midpt.x()*1000 << "," << midpt.y()*1000 << ")mm"
              << ", |w|=" << std::abs(wp1_ana.weights[i])
              << ", sin(pi*x/a)=" << expected_pattern
              << ", edge_y=" << edge_y_len*1000 << "mm" << std::endl;
    shown++;
  }
  std::cout << "Port 2 (first 5 edges):" << std::endl;
  for (int i = 0; i < std::min(5, (int)wp2_ana.edges.size()); ++i) {
    int e = wp2_ana.edges[i];
    std::cout << "  Edge " << e << ": x=" << res.x(e)
              << ", w=" << wp2_ana.weights[i]
              << ", PEC=" << (bc.dirichlet_edges.count(e) ? "yes" : "no") << std::endl;
  }

  // TEST: Solve WITHOUT port loading (just source term)
  std::cout << "\n=== Test without port loading ===" << std::endl;
  auto asmbl_no_Y = assemble_maxwell(mesh, p, bc, {}, -1);
  // Add just the source term, not the Y
  for (size_t a = 0; a < wp1_ana.edges.size(); ++a) {
    int ga = wp1_ana.edges[a];
    if (!bc.dirichlet_edges.count(ga)) {
      const std::complex<double> scale = 2.0 / std::sqrt(wp1_ana.mode.Z0);
      asmbl_no_Y.b(ga) += scale * wp1_ana.weights[a];
    }
  }
  auto res_no_Y = solve_linear(asmbl_no_Y.A, asmbl_no_Y.b, {});
  std::cout << "Solver converged (no Y): " << res_no_Y.converged << std::endl;

  // Check residual
  VecC residual_no_Y = asmbl_no_Y.A * res_no_Y.x - asmbl_no_Y.b;
  double res_norm_no_Y = residual_no_Y.norm();
  double b_norm_no_Y = asmbl_no_Y.b.norm();
  std::cout << "Residual |Ax-b|: " << res_norm_no_Y << std::endl;
  std::cout << "RHS |b|: " << b_norm_no_Y << std::endl;
  std::cout << "Relative residual: " << res_norm_no_Y / b_norm_no_Y << std::endl;

  // Check field at both ports (no Y)
  std::complex<double> V1_no_Y = 0, V2_no_Y = 0;
  for (size_t k = 0; k < wp1_ana.edges.size(); ++k) {
    int e = wp1_ana.edges[k];
    if (!bc.dirichlet_edges.count(e))
      V1_no_Y += std::conj(wp1_ana.weights[k]) * res_no_Y.x(e);
  }
  for (size_t k = 0; k < wp2_ana.edges.size(); ++k) {
    int e = wp2_ana.edges[k];
    if (!bc.dirichlet_edges.count(e))
      V2_no_Y += std::conj(wp2_ana.weights[k]) * res_no_Y.x(e);
  }
  std::cout << "V1 (no Y): " << V1_no_Y << std::endl;
  std::cout << "V2 (no Y): " << V2_no_Y << std::endl;
  std::cout << "V_inc: " << std::sqrt(wp1_ana.mode.Z0) << std::endl;
  std::cout << "|V2/V1|: " << std::abs(V2_no_Y / V1_no_Y) << std::endl;

  // Check if any interior edges have non-zero solution
  std::cout << "\n=== Interior edge analysis ===" << std::endl;
  std::set<int> surface_edges;
  for (const auto& tri : mesh.tris) {
    for (int e = 0; e < 3; ++e) {
      surface_edges.insert(tri.edges[e]);
    }
  }
  int interior_count = 0;
  double interior_sum_sq = 0, interior_max_sq = 0;
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    if (!surface_edges.count(i) && !bc.dirichlet_edges.count(i)) {
      interior_count++;
      double x_sq = std::norm(res.x(i));
      interior_sum_sq += x_sq;
      interior_max_sq = std::max(interior_max_sq, x_sq);
    }
  }
  std::cout << "Interior edges (not on any surface): " << interior_count << std::endl;
  std::cout << "Avg |x|^2 at interior: " << (interior_count > 0 ? interior_sum_sq/interior_count : 0) << std::endl;
  std::cout << "Max |x|^2 at interior: " << interior_max_sq << std::endl;
  std::cout << "For comparison, port 1 avg |x|^2: " << (1.03482e-05) << std::endl;

  // Check field distribution along z-axis
  std::cout << "\n=== Field distribution along z ===" << std::endl;
  // Group edges by z-coordinate of their midpoint
  std::map<int, std::pair<double, int>> z_bin_data; // bin -> (sum_x_sq, count)
  for (size_t i = 0; i < mesh.edges.size(); ++i) {
    if (bc.dirichlet_edges.count(i))
      continue;
    const auto& edge = mesh.edges[i];
    const auto& pa = mesh.nodes[mesh.nodeIndex.at(edge.n0)].xyz;
    const auto& pb = mesh.nodes[mesh.nodeIndex.at(edge.n1)].xyz;
    double z_mid = (pa.z() + pb.z()) / 2.0;
    int bin = static_cast<int>(z_mid * 100); // 1cm bins
    double x_sq = std::norm(res.x(i));
    if (z_bin_data.find(bin) == z_bin_data.end()) {
      z_bin_data[bin] = {0.0, 0};
    }
    z_bin_data[bin].first += x_sq;
    z_bin_data[bin].second += 1;
  }
  std::cout << "Z(cm) | Avg |x|^2 | Edge count" << std::endl;
  for (const auto& kv : z_bin_data) {
    double z = kv.first / 100.0;
    double avg_x_sq = kv.second.first / kv.second.second;
    std::cout << z*100 << " | " << avg_x_sq << " | " << kv.second.second << std::endl;
  }

  return 0;
}
