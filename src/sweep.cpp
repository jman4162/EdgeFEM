#include "edgefem/sweep.hpp"

#include "edgefem/edge_basis.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <chrono>
#include <iostream>
#include <unordered_set>

namespace edgefem {

namespace {
constexpr double c0 = 299792458.0; // speed of light (m/s)

// ---- Toy two-resonator (kept for backward compatibility) ----

Eigen::Matrix2d build_A(double f) {
  const double f0 = 1.0;
  Eigen::Matrix2d A0;
  A0 << 4.0, 1.0, 1.0, 3.0;
  Eigen::Matrix2d A1;
  A1 << 0.1, 0.05, 0.05, 0.2;
  return A0 + (f - f0) * A1;
}

Eigen::Vector2d solve_iterative(const Eigen::Matrix2d &A,
                                const Eigen::Vector2d &b,
                                const Eigen::Vector2d &x0, int &iters) {
  Eigen::Vector2d x = x0;
  const double omega = 0.1;
  iters = 0;
  for (; iters < 1000; ++iters) {
    Eigen::Vector2d r = b - A * x;
    if (r.norm() < 1e-8)
      break;
    x += omega * r;
  }
  return x;
}

} // namespace

SweepLog sweep_two_resonator(const std::vector<double> &freqs,
                             SweepPolicy policy) {
  SweepLog log;
  log.solve_time.reserve(freqs.size());
  log.iterations.reserve(freqs.size());
  log.s21.reserve(freqs.size());

  Eigen::Vector2d x_prev = Eigen::Vector2d::Zero();
  const Eigen::Vector2d b(1.0, 0.0);
  size_t reset_period = freqs.size() / 10 + 1;

  for (size_t i = 0; i < freqs.size(); ++i) {
    const double f = freqs[i];
    Eigen::Vector2d guess = Eigen::Vector2d::Zero();
    if (policy == SweepPolicy::Reuse && i > 0) {
      guess = x_prev;
    } else if (policy == SweepPolicy::Balanced && i > 0 &&
               (i % reset_period) != 0) {
      guess = x_prev;
    }

    auto A = build_A(f);
    int iters = 0;
    auto start = std::chrono::steady_clock::now();
    auto x = solve_iterative(A, b, guess, iters);
    auto end = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(end - start).count();

    log.solve_time.push_back(dt);
    log.iterations.push_back(iters);
    log.s21.push_back(std::complex<double>(x[1], 0.0));
    x_prev = x;
  }
  return log;
}

// ---- Real Maxwell frequency sweep ----

SpMatC KMMatrices::combine(double omega) const {
  const double k0 = omega / c0;
  const double k0_sq = k0 * k0;
  // A = K - k₀² M
  SpMatC A = K - k0_sq * M;
  return A;
}

KMMatrices assemble_maxwell_km(const Mesh &mesh, const MaxwellParams &p,
                               const BC &bc) {
  const int m = static_cast<int>(mesh.edges.size());

  using T = Eigen::Triplet<std::complex<double>>;
  std::vector<T> k_trips, m_trips;

  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }

    Eigen::Matrix<std::complex<double>, 6, 6> Kloc =
        whitney_curl_curl_matrix(X).cast<std::complex<double>>();
    Eigen::Matrix<std::complex<double>, 6, 6> Mloc =
        whitney_mass_matrix(X).cast<std::complex<double>>();

    // Get per-element material properties (static, no frequency dependence)
    std::complex<double> eps_r_elem = p.eps_r;
    std::complex<double> mu_r_elem = p.mu_r;
    auto eps_it = p.eps_r_regions.find(tet.phys);
    if (eps_it != p.eps_r_regions.end())
      eps_r_elem = eps_it->second;
    auto mu_it = p.mu_r_regions.find(tet.phys);
    if (mu_it != p.mu_r_regions.end())
      mu_r_elem = mu_it->second;

    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int si = tet.edge_orient[i];
      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int sj = tet.edge_orient[j];
        double sign = static_cast<double>(si * sj);

        // K contribution: (1/μ_r) * curl-curl
        std::complex<double> k_val = (Kloc(i, j) / mu_r_elem) * sign;
        k_trips.emplace_back(gi, gj, k_val);

        // M contribution: ε_r * mass
        std::complex<double> m_val = (eps_r_elem * Mloc(i, j)) * sign;
        m_trips.emplace_back(gi, gj, m_val);
      }
    }
  }

  KMMatrices km;
  km.K.resize(m, m);
  km.M.resize(m, m);
  km.K.setFromTriplets(k_trips.begin(), k_trips.end());
  km.M.setFromTriplets(m_trips.begin(), m_trips.end());
  km.K.makeCompressed();
  km.M.makeCompressed();

  // Apply Dirichlet BC to K and M: zero rows/cols, set K diagonal to 1
  for (int e : bc.dirichlet_edges) {
    // Zero column
    for (SpMatC::InnerIterator it(km.K, e); it; ++it)
      it.valueRef() = 0.0;
    for (SpMatC::InnerIterator it(km.M, e); it; ++it)
      it.valueRef() = 0.0;
  }
  // Zero rows
  for (int col = 0; col < m; ++col) {
    for (SpMatC::InnerIterator it(km.K, col); it; ++it) {
      if (bc.dirichlet_edges.count(it.row()))
        it.valueRef() = 0.0;
    }
    for (SpMatC::InnerIterator it(km.M, col); it; ++it) {
      if (bc.dirichlet_edges.count(it.row()))
        it.valueRef() = 0.0;
    }
  }
  // Set K diagonal to 1 for Dirichlet DOFs
  for (int e : bc.dirichlet_edges) {
    km.K.coeffRef(e, e) = 1.0;
    // M diagonal stays 0 so k₀²*0 = 0 and A(e,e) = 1
  }

  return km;
}

SweepResult frequency_sweep(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports,
                            const std::vector<double> &frequencies,
                            const SolveOptions &opts) {
  SweepResult result;
  result.frequencies = frequencies;
  result.S_matrices.reserve(frequencies.size());

  const int num_ports = static_cast<int>(ports.size());
  if (num_ports == 0 || frequencies.empty()) {
    return result;
  }

  // Check if we can use K/M separation (no PML, no dispersive materials)
  bool has_pml = !p.pml_regions.empty() || !p.pml_tensor_regions.empty();
  bool has_dispersive = !p.eps_models.empty() || !p.mu_models.empty();
  bool can_separate = !has_pml && !has_dispersive;

  if (can_separate) {
    // Fast path: pre-assemble K and M once
    if (opts.verbose) {
      std::cerr << "FrequencySweep: Using K/M separation ("
                << frequencies.size() << " points, " << num_ports << " ports)\n";
    }

    KMMatrices km = assemble_maxwell_km(mesh, p, bc);

    // Pre-compute port loading (frequency-independent for lumped ports)
    // Port admittance: Y_ij = w_i * conj(w_j) / Z0
    using T = Eigen::Triplet<std::complex<double>>;
    std::vector<T> port_trips;
    for (const auto &port : ports) {
      if (port.weights.size() != static_cast<int>(port.edges.size()))
        continue;
      if (port.mode.Z0 == std::complex<double>(0.0))
        continue;
      for (size_t a = 0; a < port.edges.size(); ++a) {
        int ga = port.edges[a];
        if (bc.dirichlet_edges.count(ga))
          continue;
        for (size_t b = 0; b < port.edges.size(); ++b) {
          int gb = port.edges[b];
          if (bc.dirichlet_edges.count(gb))
            continue;
          port_trips.emplace_back(
              ga, gb,
              port.weights[a] * std::conj(port.weights[b]) / port.mode.Z0);
        }
      }
    }

    const int m_size = static_cast<int>(mesh.edges.size());
    SpMatC Y_port(m_size, m_size);
    Y_port.setFromTriplets(port_trips.begin(), port_trips.end());
    Y_port.makeCompressed();

    // ABC edge set (frequency-independent geometry)
    std::unordered_set<int> abc_edges;
    if (p.use_abc) {
      for (const auto &tri : mesh.tris) {
        if (!p.abc_surface_tags.empty() && !p.abc_surface_tags.count(tri.phys))
          continue;
        for (int e = 0; e < 3; ++e) {
          int ge = tri.edges[e];
          if (!bc.dirichlet_edges.count(ge))
            abc_edges.insert(ge);
        }
      }
    }

    // Symbolic factorization (reused across frequencies for direct solver)
    Eigen::SparseLU<SpMatC> lu_solver;
    bool symbolic_done = false;

    for (size_t fi = 0; fi < frequencies.size(); ++fi) {
      double freq = frequencies[fi];
      double omega = 2.0 * M_PI * freq;
      double k0 = omega / c0;

      // Combine A = K - k₀²M + Y_port + ABC
      SpMatC A = km.combine(omega) + Y_port;

      // Add ABC: jk₀ on diagonal of ABC edges
      if (p.use_abc) {
        std::complex<double> abc_coeff(0.0, k0);
        for (int e : abc_edges) {
          A.coeffRef(e, e) += abc_coeff;
        }
      }

      // Solve for each active port
      Eigen::MatrixXcd S(num_ports, num_ports);

      for (int active = 0; active < num_ports; ++active) {
        // Build RHS: b = 2w/sqrt(Z0) for active port
        VecC b_rhs = VecC::Zero(m_size);
        const auto &aport = ports[active];
        std::complex<double> scale = 2.0 / std::sqrt(aport.mode.Z0);
        for (size_t a = 0; a < aport.edges.size(); ++a) {
          int ga = aport.edges[a];
          if (!bc.dirichlet_edges.count(ga))
            b_rhs(ga) += scale * aport.weights[a];
        }

        // Solve
        SolveResult sol;
        if (opts.use_direct) {
          if (!symbolic_done) {
            lu_solver.analyzePattern(A);
            symbolic_done = true;
          }
          lu_solver.factorize(A);
          if (lu_solver.info() != Eigen::Success) {
            sol.converged = false;
            sol.method = "SparseLU";
            sol.error_message = "Factorization failed";
          } else {
            sol.x = lu_solver.solve(b_rhs);
            sol.converged = (lu_solver.info() == Eigen::Success);
            sol.method = "SparseLU";
            sol.iters = 1;
            sol.residual = (A * sol.x - b_rhs).norm() / b_rhs.norm();
          }
        } else {
          sol = solve_linear(A, b_rhs, opts);
        }

        if (!sol.converged) {
          for (int j = 0; j < num_ports; ++j) {
            S(j, active) = std::complex<double>(
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN());
          }
          continue;
        }

        // Extract S-parameters
        std::complex<double> V_inc = std::sqrt(aport.mode.Z0);
        for (int j = 0; j < num_ports; ++j) {
          std::complex<double> Vj = 0.0;
          for (size_t k = 0; k < ports[j].edges.size(); ++k) {
            int edge_idx = ports[j].edges[k];
            if (!bc.dirichlet_edges.count(edge_idx))
              Vj += std::conj(ports[j].weights[k]) * sol.x(edge_idx);
          }
          if (j == active) {
            S(j, active) = (Vj - V_inc) / V_inc;
          } else {
            S(j, active) = Vj / V_inc;
          }
        }
      }

      result.S_matrices.push_back(S);
    }
  } else {
    // Slow path: full assembly at each frequency (dispersive/PML cases)
    if (opts.verbose) {
      std::cerr << "FrequencySweep: Using per-frequency assembly ("
                << frequencies.size()
                << " points, dispersive=" << has_dispersive
                << ", PML=" << has_pml << ")\n";
    }

    for (size_t fi = 0; fi < frequencies.size(); ++fi) {
      MaxwellParams p_freq = p;
      p_freq.omega = 2.0 * M_PI * frequencies[fi];

      Eigen::MatrixXcd S =
          calculate_sparams(mesh, p_freq, bc, ports, opts);
      result.S_matrices.push_back(S);
    }
  }

  return result;
}

} // namespace edgefem
