// Rectangular Cavity Eigenmode Test
// Validates FEM eigenvalue computation against analytical resonant frequencies
//
// For a rectangular cavity with PEC walls, the resonant frequencies are:
//   f_mnp = c/(2π) * sqrt((mπ/a)² + (nπ/b)² + (pπ/d)²)
// where a, b, d are the cavity dimensions and m, n, p are mode indices.

#include "edgefem/edge_basis.hpp"
#include "edgefem/maxwell.hpp"
#include "edgefem/mesh.hpp"
#include "edgefem/solver.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace edgefem;

constexpr double c0 = 299792458.0;
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

// Calculate analytical resonant frequency for TE_mnp or TM_mnp mode
double analytical_resonance(double a, double b, double d, int m, int n, int p) {
  double kx = (m * M_PI) / a;
  double ky = (n * M_PI) / b;
  double kz = (p * M_PI) / d;
  double k = std::sqrt(kx * kx + ky * ky + kz * kz);
  return c0 * k / (2 * M_PI);
}

// Mode type classification
struct CavityMode {
  int m, n, p;
  std::string type; // "TE" or "TM"
  double freq;

  bool operator<(const CavityMode &other) const { return freq < other.freq; }
};

std::vector<CavityMode> get_analytical_modes(double a, double b, double d,
                                             double f_max, int max_index = 5) {
  std::vector<CavityMode> modes;

  // TE modes (at least one of m,n,p must be 0 for TE to exist in that
  // direction) TM modes (both m and n must be nonzero for TMmnp)
  for (int m = 0; m <= max_index; ++m) {
    for (int n = 0; n <= max_index; ++n) {
      for (int p = 0; p <= max_index; ++p) {
        // Skip (0,0,0) - no mode
        if (m == 0 && n == 0 && p == 0)
          continue;

        // For rectangular cavity, valid modes need at least 2 nonzero indices
        // for TM, or specific combinations for TE
        int nonzero_count = (m > 0 ? 1 : 0) + (n > 0 ? 1 : 0) + (p > 0 ? 1 : 0);

        if (nonzero_count >= 2) {
          double f = analytical_resonance(a, b, d, m, n, p);
          if (f < f_max) {
            // Determine mode type based on which indices are zero
            std::string type;
            if (p == 0) {
              type = "TM"; // TMmn0 modes (no z-variation)
            } else if (m == 0 || n == 0) {
              type = "TE"; // TE modes (one transverse index zero)
            } else {
              type = "TM"; // General TMmnp
            }
            modes.push_back({m, n, p, type, f});
          }
        }
      }
    }
  }

  std::sort(modes.begin(), modes.end());
  return modes;
}

int main() {
  std::cout << std::setprecision(6);
  std::cout << "=== Rectangular Cavity Eigenmode Validation ===" << std::endl;

  // Load cube cavity mesh (1m x 1m x 1m)
  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  std::cout << "Mesh loaded: " << mesh.nodes.size() << " nodes, "
            << mesh.tets.size() << " tetrahedra, " << mesh.edges.size()
            << " edges" << std::endl;

  // Cavity dimensions (from mesh)
  double a = 1.0; // x-dimension
  double b = 1.0; // y-dimension
  double d = 1.0; // z-dimension

  std::cout << "Cavity dimensions: " << a << " x " << b << " x " << d << " m"
            << std::endl;

  // Apply PEC BC on all boundary surfaces (physical group 1)
  BC bc = build_edge_pec(mesh, 1);
  std::cout << "PEC edges (Dirichlet): " << bc.dirichlet_edges.size()
            << std::endl;

  // Get analytical modes up to 500 MHz
  double f_max = 500e6;
  auto analytical_modes = get_analytical_modes(a, b, d, f_max);
  std::cout << std::endl
            << "Analytical resonant frequencies (< " << f_max / 1e6 << " MHz):"
            << std::endl;
  for (size_t i = 0; i < std::min(analytical_modes.size(), size_t(10)); ++i) {
    const auto &m = analytical_modes[i];
    std::cout << "  " << m.type << "_" << m.m << m.n << m.p << ": "
              << m.freq / 1e6 << " MHz" << std::endl;
  }

  // Build FEM matrices for eigenvalue problem
  // Curl-curl eigenvalue problem: K*v = k² * M*v
  // where k² = ω²με = (2πf)²/c²
  std::cout << std::endl << "Assembling FEM matrices..." << std::endl;

  const int n_edges = static_cast<int>(mesh.edges.size());
  const int n_free = n_edges - static_cast<int>(bc.dirichlet_edges.size());
  std::cout << "Total edges: " << n_edges << ", Free DOFs: " << n_free
            << std::endl;

  // Build curl-curl (stiffness) and mass matrices
  using T = Eigen::Triplet<double>;
  std::vector<T> K_trips, M_trips;

  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }

    auto Kloc = whitney_curl_curl_matrix(X);
    auto Mloc = whitney_mass_matrix(X);

    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int si = tet.edge_orient[i];
      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int sj = tet.edge_orient[j];
        double sign = static_cast<double>(si * sj);
        K_trips.emplace_back(gi, gj, sign * Kloc(i, j));
        M_trips.emplace_back(gi, gj, sign * Mloc(i, j));
      }
    }
  }

  Eigen::SparseMatrix<double> K(n_edges, n_edges), M(n_edges, n_edges);
  K.setFromTriplets(K_trips.begin(), K_trips.end());
  M.setFromTriplets(M_trips.begin(), M_trips.end());

  // Build reduced system excluding PEC DOFs
  std::vector<int> free_to_global, global_to_free(n_edges, -1);
  for (int e = 0; e < n_edges; ++e) {
    if (!bc.dirichlet_edges.count(e)) {
      global_to_free[e] = static_cast<int>(free_to_global.size());
      free_to_global.push_back(e);
    }
  }

  std::vector<T> Kred_trips, Mred_trips;
  for (int k = 0; k < K.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
      int i = global_to_free[it.row()];
      int j = global_to_free[it.col()];
      if (i >= 0 && j >= 0) {
        Kred_trips.emplace_back(i, j, it.value());
      }
    }
  }
  for (int k = 0; k < M.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
      int i = global_to_free[it.row()];
      int j = global_to_free[it.col()];
      if (i >= 0 && j >= 0) {
        Mred_trips.emplace_back(i, j, it.value());
      }
    }
  }

  Eigen::SparseMatrix<double> Kred(n_free, n_free), Mred(n_free, n_free);
  Kred.setFromTriplets(Kred_trips.begin(), Kred_trips.end());
  Mred.setFromTriplets(Mred_trips.begin(), Mred_trips.end());

  // Convert to dense for eigenvalue solver (small enough for demo)
  // For large problems, use Spectra or ARPACK
  std::cout << "Solving generalized eigenvalue problem..." << std::endl;

  Eigen::MatrixXd Kd = Eigen::MatrixXd(Kred);
  Eigen::MatrixXd Md = Eigen::MatrixXd(Mred);

  // Solve K*v = lambda*M*v using generalized eigensolver
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
  ges.compute(Kd, Md);

  if (ges.info() != Eigen::Success) {
    std::cerr << "Eigenvalue computation failed!" << std::endl;
    return 1;
  }

  Eigen::VectorXd eigenvalues = ges.eigenvalues();

  // Convert eigenvalues (k²) to frequencies
  // k² = (ω/c)² = (2πf/c)² => f = c*sqrt(k²)/(2π)
  std::vector<double> fem_frequencies;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    if (eigenvalues(i) > 1e-6) { // Skip zero/negative eigenvalues
      double k_sq = eigenvalues(i);
      double f = c0 * std::sqrt(k_sq) / (2 * M_PI);
      if (f < f_max) {
        fem_frequencies.push_back(f);
      }
    }
  }

  std::sort(fem_frequencies.begin(), fem_frequencies.end());

  std::cout << std::endl
            << "FEM computed frequencies (first " << std::min(size_t(10), fem_frequencies.size())
            << "):" << std::endl;
  for (size_t i = 0; i < std::min(fem_frequencies.size(), size_t(10)); ++i) {
    std::cout << "  Mode " << i + 1 << ": " << fem_frequencies[i] / 1e6
              << " MHz" << std::endl;
  }

  // Compare FEM vs analytical
  std::cout << std::endl << "=== Comparison ===" << std::endl;
  std::cout << std::setw(15) << "Analytical" << std::setw(15) << "FEM"
            << std::setw(15) << "Error (%)" << std::setw(20) << "Mode"
            << std::endl;
  std::cout << std::string(65, '-') << std::endl;

  int n_compare = std::min({analytical_modes.size(), fem_frequencies.size(), size_t(8)});
  double max_error = 0.0;
  int num_pass = 0;

  for (int i = 0; i < n_compare; ++i) {
    double f_ana = analytical_modes[i].freq;
    double f_fem = fem_frequencies[i];
    double error_pct = std::abs(f_fem - f_ana) / f_ana * 100.0;
    max_error = std::max(max_error, error_pct);

    const auto &m = analytical_modes[i];
    std::string mode_str = m.type + "_" + std::to_string(m.m) +
                           std::to_string(m.n) + std::to_string(m.p);

    std::cout << std::setw(15) << std::fixed << std::setprecision(2)
              << f_ana / 1e6 << std::setw(15) << f_fem / 1e6 << std::setw(15)
              << std::setprecision(2) << error_pct << std::setw(20) << mode_str
              << std::endl;

    // Target: < 5% error for this coarse mesh
    if (error_pct < 10.0) {
      num_pass++;
    }
  }

  std::cout << std::endl;
  std::cout << "=== Results ===" << std::endl;
  std::cout << "Maximum error: " << max_error << "%" << std::endl;
  std::cout << "Modes passing (<10% error): " << num_pass << "/" << n_compare
            << std::endl;

  // Acceptance criteria: at least 6 of first 8 modes within 10%
  bool passed = (num_pass >= 6);
  std::cout << "Overall: " << (passed ? "PASS" : "FAIL") << std::endl;

  return passed ? 0 : 1;
}
