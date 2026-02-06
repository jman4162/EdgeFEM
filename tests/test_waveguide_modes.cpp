// Test if FEM eigenvalues match expected waveguide mode cutoffs
// The generalized eigenvalue problem K*x = λ*M*x gives λ = kc² (cutoff wavenumber squared)
// For WR-90: TE10 has kc = π/a = 137.4, so kc² = 18886

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "edgefem/maxwell.hpp"

using namespace edgefem;

constexpr double c0 = 299792458.0;

int main() {
  std::cout << std::setprecision(6);

  Mesh mesh = load_gmsh_v2("examples/rect_waveguide.msh");
  BC bc = build_edge_pec(mesh, 1);

  // Expected waveguide mode cutoffs for WR-90 (a=22.86mm, b=10.16mm)
  double a = 0.02286, b = 0.01016, L = 0.05;

  std::cout << "=== Expected Waveguide Mode Cutoffs ===" << std::endl;
  std::cout << "WR-90: a=" << a*1000 << "mm, b=" << b*1000 << "mm, L=" << L*1000 << "mm" << std::endl;

  // TE_mn modes: kc = π*sqrt((m/a)² + (n/b)²)
  // TM_mn modes: same formula (but TM01/TM10 don't exist for rect WG)
  std::vector<std::tuple<int,int,double,std::string>> modes;
  for (int m = 0; m <= 3; ++m) {
    for (int n = 0; n <= 3; ++n) {
      if (m == 0 && n == 0) continue;  // No TE00 or TM00
      double kc = M_PI * std::sqrt(std::pow(m/a, 2) + std::pow(n/b, 2));
      double fc = kc * c0 / (2 * M_PI);
      std::string name = "TE" + std::to_string(m) + std::to_string(n);
      modes.push_back({m, n, fc, name});

      // TM modes (m,n both nonzero)
      if (m > 0 && n > 0) {
        std::string tm_name = "TM" + std::to_string(m) + std::to_string(n);
        modes.push_back({m, n, fc, tm_name});
      }
    }
  }

  std::sort(modes.begin(), modes.end(),
            [](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });

  std::cout << "\nFirst 10 modes by cutoff frequency:" << std::endl;
  for (int i = 0; i < std::min(10, (int)modes.size()); ++i) {
    auto& [m, n, fc, name] = modes[i];
    double kc = 2*M_PI*fc/c0;
    std::cout << "  " << name << ": fc = " << fc/1e9 << " GHz, kc² = " << kc*kc << std::endl;
  }

  // Also compute cavity resonances (for waveguide with PEC at both ends)
  std::cout << "\n=== Expected Cavity Resonances ===" << std::endl;
  std::vector<std::tuple<int,int,int,double,std::string>> cavity_modes;
  for (int m = 0; m <= 3; ++m) {
    for (int n = 0; n <= 3; ++n) {
      for (int p = 0; p <= 5; ++p) {
        if (m == 0 && n == 0) continue;
        double f = (c0/2) * std::sqrt(std::pow(m/a, 2) + std::pow(n/b, 2) + std::pow(p/L, 2));
        std::string name = "TE" + std::to_string(m) + std::to_string(n) + std::to_string(p);
        cavity_modes.push_back({m, n, p, f, name});
      }
    }
  }

  std::sort(cavity_modes.begin(), cavity_modes.end(),
            [](const auto& a, const auto& b) { return std::get<3>(a) < std::get<3>(b); });

  std::cout << "First 10 cavity modes:" << std::endl;
  for (int i = 0; i < std::min(10, (int)cavity_modes.size()); ++i) {
    auto& [m, n, p, f, name] = cavity_modes[i];
    std::cout << "  " << name << ": f = " << f/1e9 << " GHz" << std::endl;
  }

  // Assemble K and M matrices separately
  std::cout << "\n=== FEM Matrix Assembly ===" << std::endl;

  // Very low frequency to get K
  MaxwellParams p_low;
  p_low.omega = 1e-10;
  auto asmbl_low = assemble_maxwell(mesh, p_low, bc, {}, -1);

  // Moderate frequency to get K - k0²M, then extract M
  double f_test = 1e9;  // 1 GHz
  double k0_test = 2*M_PI*f_test/c0;
  MaxwellParams p_test;
  p_test.omega = 2*M_PI*f_test;
  auto asmbl_test = assemble_maxwell(mesh, p_test, bc, {}, -1);

  int n = asmbl_low.A.rows();
  int n_free = n - bc.dirichlet_edges.size();

  // Extract free-DOF submatrices
  std::vector<int> free_to_orig;
  for (int i = 0; i < n; ++i) {
    if (!bc.dirichlet_edges.count(i)) {
      free_to_orig.push_back(i);
    }
  }

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

  std::cout << "Free DOFs: " << n_free << std::endl;
  std::cout << "K trace: " << K_free.trace() << std::endl;
  std::cout << "M trace: " << M_free.trace() << std::endl;

  // Solve generalized eigenvalue problem K*x = λ*M*x
  std::cout << "\nComputing generalized eigenvalues (K*x = λ*M*x)..." << std::endl;

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(K_free, M_free);

  if (ges.info() != Eigen::Success) {
    std::cout << "Eigenvalue computation failed!" << std::endl;
    return 1;
  }

  auto eigenvalues = ges.eigenvalues();

  // Sort and filter positive eigenvalues (these are kc² values)
  std::vector<double> kc_squared;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    if (eigenvalues(i) > 100) {  // Filter out near-zero (null space)
      kc_squared.push_back(eigenvalues(i));
    }
  }
  std::sort(kc_squared.begin(), kc_squared.end());

  std::cout << "\nFirst 20 FEM eigenvalues (kc²):" << std::endl;
  std::cout << "FEM kc²       FEM fc(GHz)   Expected" << std::endl;

  // Compare with expected
  double kc_te10 = M_PI / a;
  double kc_te10_sq = kc_te10 * kc_te10;
  double kc_te20 = 2 * M_PI / a;
  double kc_te20_sq = kc_te20 * kc_te20;
  double kc_te01 = M_PI / b;
  double kc_te01_sq = kc_te01 * kc_te01;

  for (int i = 0; i < std::min(20, (int)kc_squared.size()); ++i) {
    double kc_sq = kc_squared[i];
    double kc = std::sqrt(kc_sq);
    double fc = kc * c0 / (2 * M_PI);

    std::string expected = "";
    if (std::abs(kc_sq - kc_te10_sq) / kc_te10_sq < 0.1) expected = "~TE10";
    else if (std::abs(kc_sq - kc_te20_sq) / kc_te20_sq < 0.1) expected = "~TE20";
    else if (std::abs(kc_sq - kc_te01_sq) / kc_te01_sq < 0.1) expected = "~TE01";

    std::cout << std::setw(12) << kc_sq << "  " << std::setw(10) << fc/1e9 << "  " << expected << std::endl;
  }

  std::cout << "\nExpected TE10: kc² = " << kc_te10_sq << ", fc = " << kc_te10*c0/(2*M_PI)/1e9 << " GHz" << std::endl;
  std::cout << "Expected TE20: kc² = " << kc_te20_sq << ", fc = " << kc_te20*c0/(2*M_PI)/1e9 << " GHz" << std::endl;
  std::cout << "Expected TE01: kc² = " << kc_te01_sq << ", fc = " << kc_te01*c0/(2*M_PI)/1e9 << " GHz" << std::endl;

  // Check if TE10 mode is captured
  double min_diff = 1e100;
  for (double kc_sq : kc_squared) {
    min_diff = std::min(min_diff, std::abs(kc_sq - kc_te10_sq));
  }
  std::cout << "\nClosest FEM eigenvalue to TE10 kc²: diff = " << min_diff
            << " (" << min_diff/kc_te10_sq*100 << "% error)" << std::endl;

  if (min_diff / kc_te10_sq > 0.1) {
    std::cout << "WARNING: FEM may not be capturing TE10 mode accurately!" << std::endl;
  } else {
    std::cout << "TE10 mode is captured with acceptable accuracy." << std::endl;
  }

  return 0;
}
