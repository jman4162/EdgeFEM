#include "edgefem/sweep.hpp"

#include "edgefem/edge_basis.hpp"
#include "edgefem/free_dofs.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

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

namespace {

// Collects the field solutions a SweepFieldRequest asks for while the sweep
// runs. For min_reflection only the best fields so far are held.
class FieldKeeper {
public:
  FieldKeeper(const SweepFieldRequest &req, std::size_t num_freqs)
      : min_reflection_(req.min_reflection) {
    for (std::size_t i : req.indices) {
      if (i >= num_freqs)
        throw std::invalid_argument(
            "frequency_sweep: kept field index " + std::to_string(i) +
            " is out of range for " + std::to_string(num_freqs) +
            " frequencies");
      indices_.insert(i);
    }
  }

  // Offer the solutions at frequency index fi; S is its S-matrix. Takes
  // the vectors from sols only when they are kept.
  void offer(std::size_t fi, const Eigen::MatrixXcd &S,
             std::vector<VecC> &sols) {
    const bool listed = indices_.count(fi) > 0;
    bool best = false;
    if (min_reflection_ && S.size() > 0) {
      const double r = std::abs(S(0, 0));
      if (!std::isnan(r) && (!have_best_ || r < best_r_)) {
        best = true;
        have_best_ = true;
        best_r_ = r;
        best_index_ = fi;
        best_x_.clear(); // the previous best is no longer needed
      }
    }
    if (listed)
      kept_[fi] = std::move(sols);
    else if (best)
      best_x_ = std::move(sols);
  }

  /// Frequency index with the smallest |S(0,0)|, if min_reflection was
  /// requested and some frequency had a finite S(0,0).
  std::optional<std::size_t> min_reflection_index() const {
    if (have_best_)
      return best_index_;
    return std::nullopt;
  }

  std::vector<SweepField> take() {
    if (have_best_ && !kept_.count(best_index_))
      kept_[best_index_] = std::move(best_x_);
    std::vector<SweepField> out;
    out.reserve(kept_.size());
    for (auto &kv : kept_)
      out.push_back(SweepField{kv.first, std::move(kv.second)});
    return out;
  }

private:
  bool min_reflection_;
  std::set<std::size_t> indices_;
  std::map<std::size_t, std::vector<VecC>> kept_;
  bool have_best_ = false;
  double best_r_ = 0.0;
  std::size_t best_index_ = 0;
  std::vector<VecC> best_x_;
};

} // namespace

SweepResult frequency_sweep(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports,
                            const std::vector<double> &frequencies,
                            const SolveOptions &opts,
                            const SweepFieldRequest &keep) {
  SweepResult result;
  result.frequencies = frequencies;
  result.S_matrices.reserve(frequencies.size());
  FieldKeeper keeper(keep, frequencies.size());

  const int num_ports = static_cast<int>(ports.size());
  if (num_ports == 0 || frequencies.empty()) {
    return result;
  }
  warn_modal_sparams(ports);

  // Check if we can use K/M separation (no PML, no dispersive materials)
  bool has_pml = !p.pml_regions.empty() || !p.pml_tensor_regions.empty();
  bool has_dispersive = !p.eps_models.empty() || !p.mu_models.empty();
  bool can_separate = !has_pml && !has_dispersive;

  if (can_separate) {
    // Fast path: pre-assemble K and M once
    if (opts.verbose) {
      std::cerr << "FrequencySweep: Using K/M separation ("
                << frequencies.size() << " points, " << num_ports
                << " ports)\n";
    }

    // K, M and the port/ABC terms are restricted to the free (non-PEC)
    // edges once; every frequency factors the reduced matrix and the
    // solutions are expanded back to one entry per mesh edge.
    const int m_size = static_cast<int>(mesh.edges.size());
    const FreeDofs dofs(m_size, bc);
    KMMatrices km;
    {
      const KMMatrices km_full = assemble_maxwell_km(mesh, p, bc);
      km.K = dofs.restrict_matrix(km_full.K);
      km.M = dofs.restrict_matrix(km_full.M);
    }

    // Frequency-independent part of the j k0 terms: the gap-port load
    // (eta0 / Z0) v v^H of the lumped ports (derivation in
    // assemble_maxwell.cpp) plus the ABC surface mass matrix. Modal ports
    // add the frequency-independent 1.0.0 load w w^H / Z0 directly.
    SpMatC loss = port_load_matrix(ports, bc, m_size);
    if (p.use_abc)
      loss += abc_surface_matrix(mesh, p, bc);
    const SpMatC loss_free = dofs.restrict_matrix(loss);
    const SpMatC modal_load_free =
        dofs.restrict_matrix(modal_port_load_matrix(ports, bc, m_size));

    // Symbolic analysis on the first frequency, reused for the rest.
    std::unique_ptr<DirectLU> lu;
    if (opts.use_direct)
      lu = std::make_unique<DirectLU>(opts.direct_backend);

    for (size_t fi = 0; fi < frequencies.size(); ++fi) {
      double freq = frequencies[fi];
      double omega = 2.0 * M_PI * freq;
      double k0 = omega / c0;

      // A = K - k0^2 M + j k0 (L_port + ABC) + Y_modal, free edges only
      SpMatC A = km.combine(omega) + std::complex<double>(0.0, k0) * loss_free +
                 modal_load_free;
      A.makeCompressed();

      // One matrix for all ports; port `active` only changes the RHS.
      std::vector<VecC> rhs;
      rhs.reserve(ports.size());
      for (const auto &port : ports)
        rhs.push_back(
            dofs.restrict_vector(port_excitation(port, bc, m_size, k0)));

      // Direct: numeric factorization once per frequency, one solve per port.
      std::vector<SolveResult> sols = opts.use_direct
                                          ? lu->factorize_and_solve(A, rhs)
                                          : solve_linear_multi(A, rhs, opts);

      // Full-length fields; a failed solve gives a NaN column and no field.
      Eigen::MatrixXcd S(num_ports, num_ports);
      std::vector<VecC> xs(sols.size());
      for (int active = 0; active < num_ports; ++active) {
        if (!check_solver_convergence(
                sols[active], "frequency_sweep (f = " + std::to_string(freq) +
                                  " Hz, active port " + std::to_string(active) +
                                  ")")) {
          for (int j = 0; j < num_ports; ++j) {
            S(j, active) =
                std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN());
          }
          continue;
        }
        xs[active] = dofs.expand(sols[active].x);
        fill_sparam_column(S, active, ports, bc, xs[active]);
      }

      keeper.offer(fi, S, xs);
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

      SParamSolution sol =
          calculate_sparams_with_fields(mesh, p_freq, bc, ports, opts);
      keeper.offer(fi, sol.S, sol.x);
      result.S_matrices.push_back(sol.S);
    }
  }

  result.min_reflection_index = keeper.min_reflection_index();
  result.fields = keeper.take();
  return result;
}

} // namespace edgefem
