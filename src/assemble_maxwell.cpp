#include "edgefem/edge_basis.hpp"
#include "edgefem/free_dofs.hpp"
#include "edgefem/maxwell.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "edgefem/solver.hpp"

namespace edgefem {

namespace {
/// Check solver convergence and warn to stderr if it failed.
/// Returns true if converged, false otherwise.
bool check_solver_convergence(const SolveResult &res,
                              const std::string &context) {
  if (!res.converged) {
    std::cerr << "WARNING: Solver did not converge in " << context
              << " (method=" << res.method << ", iters=" << res.iters
              << ", residual=" << res.residual;
    if (!res.error_message.empty()) {
      std::cerr << ", error: " << res.error_message;
    }
    std::cerr << "). S-parameters may be unreliable." << std::endl;
    return false;
  }
  return true;
}
} // namespace

namespace {
constexpr double c0 = 299792458.0;        // speed of light in vacuum (m/s)
constexpr double mu0 = 4.0 * M_PI * 1e-7; // permeability of free space (H/m)
constexpr double eps0 =
    1.0 / (mu0 * c0 * c0);        // permittivity of free space (F/m)
constexpr double eta0 = mu0 * c0; // free-space impedance ≈ 377 ohms
} // namespace

// ===========================================================================
// Gap (lumped) port formulation
// ===========================================================================
// Time convention exp(+j omega t), as in the rest of the code: the ABC adds
// +j k0, a lossy dielectric has Im(eps_r) < 0, and the eigenmode path
// expects phase(S21) = -beta L.
//
// Field equation, scaled by mu0 like the assembled matrix:
//   curl(mu_r^-1 curl E) - k0^2 eps_r E = -j omega mu0 J = -j k0 eta0 J.
// Galerkin with edge functions N_i gives A0 x = b with
//   A0 = K / mu_r - k0^2 eps_r M (+ ABC),  b_i = -j k0 eta0 int N_i . J dV.
//
// Port geometry: a strip S of width W across which the field points along
// the unit vector e (gap length h along e). Let w_i = int_S N_i . e dS (the
// surface-integral weights of build_lumped_port). For E uniform on the
// strip, w^T x = int_S E . e dS = (V / h) h W = V W, so the gap voltage is
//   V = w^T x / W = v^T x,   v = w / W.
// build_lumped_port stores v in WavePort::weights.
//
// Source: a sheet current I entering the + terminal flows along -e, i.e.
// J = -(I / W) e on S, and delivers P = Re(V conj(I)) / 2. Then
//   b = +j k0 eta0 (I / W) w = j k0 eta0 I v.
// Load: a resistor Z0 across the gap carries V / Z0 along +e, which adds
//   -j k0 eta0 (V / Z0) v = -(j k0 eta0 / Z0) v v^T x
// to b. Moved to the left-hand side:
//   (A0 + (j k0 eta0 / (Z0 W^2)) w w^T) x = j k0 eta0 (I / W) w.     (1)
// Without the load, Z = V / I = j k0 eta0 v^T A0^-1 v
//                         = j k0 eta0 w^T A0^-1 w / W^2,
// the readout of the 2026-10-12 investigation (scripts/e2_port_algebra.py).
//
// Passivity: A0 = R + j B with R, B real symmetric and B >= 0 (ABC, and
// k0^2 |Im eps_r| M for loss). For x = A0^-1 v, v^H x = x^H R x - j x^H B x,
// so Re Z = k0 eta0 x^H B x >= 0. The same holds for c^H Z c with any port
// combination, so the N-port Z (and S) is passive. In a lossless closed
// region B = 0 and Z is purely reactive: |S11| = 1.
//
// Complex weights: the code uses v v^H for the load and V = v^H x. For real
// v (every lumped port) these equal v v^T and v^T x above. With v v^H a
// constant phase on a port's weights (eigenvector wave ports carry a
// factor j) cancels, and the load stays passive for any weights.
//
// Excitation and S: a Norton current I = 2 / sqrt(Z0) in parallel with the
// load Z0 equals a Thevenin source 2 sqrt(Z0) behind Z0, i.e. an incident
// power wave a = 1 (V+ = sqrt(Z0)). With the other ports terminated in
// their Z0 and not driven,
//   S(j, i) = V_j / sqrt(Z0_j) - delta_ij.
// For one port, V = Z Z0 I / (Z + Z0) gives S11 = (Z - Z0) / (Z + Z0).
// All N ports share the matrix of (1) with every load present; port i only
// changes the right-hand side, so one factorization serves N solves. The S
// definition assumes real Z0 (a complex Z0 enters as 1/Z0 and sqrt(Z0)).
// The load is frequency dependent only through the factor j k0:
// port_load_matrix returns L = (eta0 / Z0) v v^H and callers add j k0 L.

namespace {
bool port_has_load(const WavePort &port) {
  return port.weights.size() == static_cast<int>(port.edges.size()) &&
         port.mode.Z0 != std::complex<double>(0.0);
}
} // namespace

SpMatC port_load_matrix(const std::vector<WavePort> &ports, const BC &bc,
                        int num_edges) {
  using T = Eigen::Triplet<std::complex<double>>;
  std::vector<T> trips;
  for (const auto &port : ports) {
    if (!port_has_load(port))
      continue;
    const std::complex<double> scale = eta0 / port.mode.Z0;
    for (size_t a = 0; a < port.edges.size(); ++a) {
      const int ga = port.edges[a];
      if (bc.dirichlet_edges.count(ga))
        continue;
      for (size_t b = 0; b < port.edges.size(); ++b) {
        const int gb = port.edges[b];
        if (bc.dirichlet_edges.count(gb))
          continue;
        trips.emplace_back(
            ga, gb, scale * port.weights[a] * std::conj(port.weights[b]));
      }
    }
  }
  SpMatC L(num_edges, num_edges);
  L.setFromTriplets(trips.begin(), trips.end());
  L.makeCompressed();
  return L;
}

SpMatC abc_surface_matrix(const Mesh &mesh, const MaxwellParams &p,
                          const BC &bc) {
  const int m = static_cast<int>(mesh.edges.size());
  if (!p.use_abc)
    return SpMatC(m, m);
  // Weak form of n x curl E + j k0 n x n x E = 0 on the boundary:
  // + j k0 int_S (n x N_i) . (n x N_j) dS = j k0 int_S N_i . N_j dS for
  // the tangential traces, the triangle Whitney mass matrix. (Before
  // 2026-10, the code added j k0 to the diagonal of each ABC edge instead.)
  SpMatC Ms =
      assemble_surface_mass(mesh, p.abc_surface_tags, bc.dirichlet_edges)
          .cast<std::complex<double>>();
  Ms.makeCompressed();
  return Ms;
}

VecC port_excitation(const WavePort &port, const BC &bc, int num_edges,
                     double k0) {
  VecC b = VecC::Zero(num_edges);
  if (!port_has_load(port))
    return b;
  const std::complex<double> current = 2.0 / std::sqrt(port.mode.Z0);
  const std::complex<double> scale =
      std::complex<double>(0.0, k0 * eta0) * current;
  for (size_t a = 0; a < port.edges.size(); ++a) {
    const int ga = port.edges[a];
    if (bc.dirichlet_edges.count(ga))
      continue;
    b(ga) += scale * port.weights[a];
  }
  return b;
}

std::complex<double> port_voltage(const WavePort &port, const BC &bc,
                                  const VecC &x) {
  std::complex<double> v = 0.0;
  if (port.weights.size() != static_cast<int>(port.edges.size()))
    return v;
  for (size_t k = 0; k < port.edges.size(); ++k) {
    const int e = port.edges[k];
    if (bc.dirichlet_edges.count(e))
      continue;
    v += std::conj(port.weights[k]) * x(e);
  }
  return v;
}

void fill_sparam_column(Eigen::MatrixXcd &S, int active,
                        const std::vector<WavePort> &ports, const BC &bc,
                        const VecC &x) {
  for (int j = 0; j < static_cast<int>(ports.size()); ++j) {
    if (!port_has_load(ports[j])) {
      S(j, active) =
          std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    const std::complex<double> b_out =
        port_voltage(ports[j], bc, x) / std::sqrt(ports[j].mode.Z0);
    S(j, active) = (j == active) ? b_out - 1.0 : b_out;
  }
}

MaxwellAssembly assemble_maxwell(const Mesh &mesh, const MaxwellParams &p,
                                 const BC &bc,
                                 const std::vector<WavePort> &ports,
                                 int active_port_idx) {
  const int m = static_cast<int>(mesh.edges.size());
  MaxwellAssembly asmbl;
  asmbl.A.resize(m, m);
  asmbl.b = VecC::Zero(m);

  using T = Eigen::Triplet<std::complex<double>>;
  std::vector<T> trips;

  struct Bounds {
    Eigen::Vector3d min =
        Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d max =
        Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
    bool initialized = false;
  };

  std::unordered_map<int, Bounds> region_bounds;
  for (const auto &kv : p.pml_tensor_regions)
    region_bounds.emplace(kv.first, Bounds{});

  if (!region_bounds.empty()) {
    for (const auto &tet : mesh.tets) {
      auto it = region_bounds.find(tet.phys);
      if (it == region_bounds.end())
        continue;
      auto &bounds = it->second;
      for (int i = 0; i < 4; ++i) {
        int idx = mesh.nodeIndex.at(tet.conn[i]);
        const Eigen::Vector3d &xyz = mesh.nodes[idx].xyz;
        if (!bounds.initialized) {
          bounds.min = xyz;
          bounds.max = xyz;
          bounds.initialized = true;
        } else {
          bounds.min = bounds.min.cwiseMin(xyz);
          bounds.max = bounds.max.cwiseMax(xyz);
        }
      }
    }
  }

  asmbl.diagnostics.clear();
  asmbl.diagnostics.reserve(p.pml_tensor_regions.size());
  for (const auto &kv : p.pml_tensor_regions) {
    PMLDiagnostic diag;
    diag.region_tag = kv.first;
    diag.sigma_max = kv.second.sigma_max;
    diag.thickness = kv.second.thickness;
    const double omega_mag = std::abs(p.omega);
    for (int axis = 0; axis < 3; ++axis) {
      const double sigma_max = kv.second.sigma_max[axis];
      const double thickness = kv.second.thickness[axis];
      if (sigma_max <= 0.0 || thickness <= 0.0 || omega_mag <= 0.0) {
        diag.reflection_est[axis] = 1.0;
        continue;
      }
      const double integral =
          sigma_max * thickness / (kv.second.grading_order + 1.0);
      const double exponent = -2.0 * integral / omega_mag;
      diag.reflection_est[axis] = std::exp(exponent);
    }
    asmbl.diagnostics.push_back(diag);
  }

  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }
    // Determine if this element lies in a PML region
    std::complex<double> stretch_scalar(1.0, 0.0);
    std::array<std::complex<double>, 3> stretch_tensor{
        std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0),
        std::complex<double>(1.0, 0.0)};

    if (p.omega != 0.0) {
      auto tensor_it = p.pml_tensor_regions.find(tet.phys);
      if (tensor_it != p.pml_tensor_regions.end()) {
        const auto &spec = tensor_it->second;
        Eigen::Vector3d sigma = Eigen::Vector3d::Zero();
        Eigen::Vector3d centroid = (X[0] + X[1] + X[2] + X[3]) / 4.0;
        auto bounds_it = region_bounds.find(tet.phys);
        Eigen::Vector3d bbox_min = Eigen::Vector3d::Zero();
        Eigen::Vector3d bbox_max = Eigen::Vector3d::Zero();
        if (bounds_it != region_bounds.end() && bounds_it->second.initialized) {
          bbox_min = bounds_it->second.min;
          bbox_max = bounds_it->second.max;
        }
        constexpr double kMinProfile = 1e-3;
        for (int axis = 0; axis < 3; ++axis) {
          const double sigma_max = spec.sigma_max[axis];
          const double thickness = spec.thickness[axis];
          if (sigma_max <= 0.0 || thickness <= 0.0) {
            sigma[axis] = 0.0;
            continue;
          }
          double dist_to_min = centroid[axis] - bbox_min[axis];
          double dist_to_max = bbox_max[axis] - centroid[axis];
          double min_dist = std::min(dist_to_min, dist_to_max);
          double xi = 0.0;
          if (min_dist < thickness) {
            xi = 1.0 - (min_dist / thickness);
          }
          if (p.enforce_pml_heuristics && xi > 0.0) {
            xi = std::max(xi, kMinProfile);
          }
          sigma[axis] = sigma_max * std::pow(xi, spec.grading_order);
          if (p.enforce_pml_heuristics && sigma_max > 0.0 &&
              sigma[axis] < sigma_max * kMinProfile) {
            sigma[axis] = sigma_max * kMinProfile;
          }
          stretch_tensor[axis] =
              std::complex<double>(1.0, sigma[axis] / p.omega);
        }
        std::complex<double> stretch_sum =
            stretch_tensor[0] + stretch_tensor[1] + stretch_tensor[2];
        stretch_scalar = stretch_sum / 3.0;
      } else if (p.pml_regions.count(tet.phys)) {
        stretch_scalar = std::complex<double>(1.0, p.pml_sigma / p.omega);
        stretch_tensor = {stretch_scalar, stretch_scalar, stretch_scalar};
      }
    }
    Eigen::Matrix<std::complex<double>, 6, 6> Kloc =
        whitney_curl_curl_matrix(X).cast<std::complex<double>>() /
        stretch_scalar;
    Eigen::Matrix<std::complex<double>, 6, 6> Mloc =
        whitney_mass_matrix(X).cast<std::complex<double>>() * stretch_scalar;
    // Curl-curl equation using normalized formulation:
    // (1/μ_r)∇×∇×E - k₀²ε_r E = source
    // where k₀ = ω/c₀ is the free-space wavenumber.
    // This produces matrix entries in the range ~10² to ~10⁴,
    // which is compatible with properly normalized port weights.
    // Matrix entry = (1/μ_r)K - k₀²ε_r M
    const double k0 = p.omega / c0;
    const double k0_sq = k0 * k0;

    // Get per-element material properties (region-based or global fallback)
    // Use frequency-aware accessor to support dispersive materials
    const std::complex<double> eps_r_elem = p.get_eps_r(tet.phys, p.omega);
    const std::complex<double> mu_r_elem = p.get_mu_r(tet.phys, p.omega);

    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int si = tet.edge_orient[i];
      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int sj = tet.edge_orient[j];
        std::complex<double> val =
            (Kloc(i, j) / mu_r_elem) - (k0_sq * eps_r_elem * Mloc(i, j));
        val *= static_cast<double>(si * sj);
        trips.emplace_back(gi, gj, val);
      }
    }
  }
  asmbl.A.setFromTriplets(trips.begin(), trips.end());

  asmbl.A.makeCompressed();

  // Gap-port loads j k0 (eta0 / Z0) v v^H on every port, and the active
  // port's source j k0 eta0 (2 / sqrt(Z0)) v (derivation above). PEC edges
  // are skipped: their field is zero by definition.
  if (!ports.empty()) {
    const double k0_port = p.omega / c0;
    const SpMatC L = port_load_matrix(ports, bc, m);
    if (L.nonZeros() > 0) {
      asmbl.A += std::complex<double>(0.0, k0_port) * L;
      asmbl.A.makeCompressed();
    }
    if (active_port_idx >= 0 &&
        active_port_idx < static_cast<int>(ports.size())) {
      asmbl.b += port_excitation(ports[active_port_idx], bc, m, k0_port);
    }
  }

  if (p.use_abc) {
    // First-order ABC: j k0 times the boundary surface mass matrix. With the
    // normalized formulation (K/mu_r - k0^2 eps_r M) the coefficient is
    // j k0 (wavenumber), not j omega.
    const double k0_abc = p.omega / c0;
    asmbl.A +=
        std::complex<double>(0.0, k0_abc) * abc_surface_matrix(mesh, p, bc);
    asmbl.A.makeCompressed();
  }

  // Port-specific ABC: Add absorbing boundary conditions on port surfaces
  // For Maxwell's equations, the first-order ABC is: n×(∇×E) + jβ E_t = 2jβ
  // E_inc In weak form this adds jβ ∫∫_S (n×N_i)·(n×N_j) dS to the stiffness
  // matrix where N_i are the edge basis functions.
  //
  // For edge elements, the surface integral ∫∫ (n×N_i)·(n×N_j) dS
  // equals the edge mass matrix on the 2D surface.
  if (p.use_port_abc && p.port_abc_type != PortABCType::None) {
    for (const auto &port : ports) {
      if (port.mode.Z0 == std::complex<double>(0.0))
        continue;

      // Compute propagation constant beta = sqrt(k0² - kc²)
      const double k0 = p.omega / c0;
      const double kc = port.mode.kc;
      const double beta_sq = k0 * k0 - kc * kc;

      if (beta_sq > 0) {
        const double beta = std::sqrt(beta_sq);
        const double Z0_real = std::real(port.mode.Z0);
        std::complex<double> abc_coeff(0.0, 0.0);

        switch (p.port_abc_type) {
        case PortABCType::Beta:
          // Standard ABC: jβ (propagation constant)
          abc_coeff = std::complex<double>(0.0, beta);
          break;
        case PortABCType::BetaNorm:
          // Dimensionless: jβ/k0
          abc_coeff = std::complex<double>(0.0, beta / k0);
          break;
        case PortABCType::ImpedanceMatch:
          // Impedance matched: jβ√(Z0/η0)
          abc_coeff =
              std::complex<double>(0.0, beta * std::sqrt(Z0_real / eta0));
          break;
        case PortABCType::ModalAdmittance:
          // Modal admittance: jωε0/Z0 = jβ/(η0·μr) for TE modes
          abc_coeff = std::complex<double>(0.0, p.omega * eps0 / Z0_real);
          break;
        case PortABCType::None:
          break;
        }

        // Add ABC contribution as diagonal term (simplified)
        // A proper implementation would use surface mass matrix integrals
        for (int e : port.edges) {
          if (!bc.dirichlet_edges.count(e)) {
            asmbl.A.coeffRef(e, e) += abc_coeff;
          }
        }
      }
    }
  }

  // Apply Dirichlet BC: zero both row and column, set diagonal to 1
  // First pass: zero column entries (column-major iteration)
  for (int e : bc.dirichlet_edges) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A,
                                                                     e);
         it; ++it) {
      it.valueRef() = std::complex<double>(0.0, 0.0);
    }
  }

  // Second pass: zero row entries by iterating over all columns
  for (int col = 0; col < m; ++col) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A,
                                                                     col);
         it; ++it) {
      if (bc.dirichlet_edges.count(it.row())) {
        it.valueRef() = std::complex<double>(0.0, 0.0);
      }
    }
  }

  // Set diagonal to 1 and RHS to 0 for Dirichlet DOFs
  for (int e : bc.dirichlet_edges) {
    asmbl.A.coeffRef(e, e) = std::complex<double>(1.0, 0.0);
    asmbl.b(e) = std::complex<double>(0.0, 0.0);
  }

  return asmbl;
}

SParamSolution
calculate_sparams_with_fields(const Mesh &mesh, const MaxwellParams &p,
                              const BC &bc, const std::vector<WavePort> &ports,
                              const SolveOptions &opts) {
  const int num_ports = static_cast<int>(ports.size());
  const int m = static_cast<int>(mesh.edges.size());
  SParamSolution out;
  out.S.resize(num_ports, num_ports);
  if (num_ports == 0)
    return out;

  // One loaded matrix for all ports; port i only changes the RHS. The
  // factored system has the Dirichlet rows and columns removed (they are
  // identity rows in asmbl.A with zero RHS, so the free-edge solution is
  // the same); solutions are expanded back to one entry per edge.
  const auto asmbl = assemble_maxwell(mesh, p, bc, ports, -1);
  const FreeDofs dofs(m, bc);
  const double k0 = p.omega / c0;
  std::vector<VecC> rhs;
  rhs.reserve(ports.size());
  for (const auto &port : ports)
    rhs.push_back(dofs.restrict_vector(port_excitation(port, bc, m, k0)));
  std::vector<SolveResult> sols =
      solve_linear_multi(dofs.restrict_matrix(asmbl.A), rhs, opts);

  out.x.reserve(ports.size());
  for (int i = 0; i < num_ports; ++i) { // active port
    if (!check_solver_convergence(sols[i], "calculate_sparams (active port " +
                                               std::to_string(i) + ")")) {
      // Fill column with NaN to signal bad results; no field.
      for (int j = 0; j < num_ports; ++j) {
        out.S(j, i) =
            std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN());
      }
      out.x.emplace_back();
      continue;
    }
    VecC x = dofs.expand(sols[i].x);
    fill_sparam_column(out.S, i, ports, bc, x);
    out.x.push_back(std::move(x));
  }
  return out;
}

Eigen::MatrixXcd calculate_sparams(const Mesh &mesh, const MaxwellParams &p,
                                   const BC &bc,
                                   const std::vector<WavePort> &ports,
                                   const SolveOptions &opts) {
  return calculate_sparams_with_fields(mesh, p, bc, ports, opts).S;
}

void normalize_port_weights(const Mesh & /*mesh*/, const MaxwellParams & /*p*/,
                            const BC & /*bc*/,
                            std::vector<WavePort> & /*ports*/,
                            const SolveOptions & /*opts*/) {
  static std::atomic<bool> warned{false};
  if (!warned.exchange(true)) {
    std::cerr << "WARNING: normalize_port_weights is deprecated and leaves the "
                 "ports unchanged. calculate_sparams and frequency_sweep use "
                 "the gap-port load j k0 eta0 / Z0 v v^H, which needs no "
                 "weight normalization. Remove the call; the function will be "
                 "removed in a later release."
              << std::endl;
  }
}

MaxwellAssembly assemble_maxwell_periodic(const Mesh &mesh,
                                          const MaxwellParams &p, const BC &bc,
                                          const PeriodicBC &pbc,
                                          const std::vector<WavePort> &ports,
                                          int active_port_idx) {
  // First, assemble the standard Maxwell system
  MaxwellAssembly asmbl = assemble_maxwell(mesh, p, bc, ports, active_port_idx);

  if (pbc.pairs.empty()) {
    return asmbl;
  }

  const std::complex<double> phi = pbc.phase_shift;
  const std::complex<double> phi_conj = std::conj(phi);

  // Build set of slave edges for quick lookup
  std::unordered_set<int> slave_edges;
  for (const auto &pair : pbc.pairs) {
    slave_edges.insert(pair.slave_edge);
  }

  // Apply periodic constraint by direct elimination:
  // For each pair, we have: E_slave = phi * orient * E_master
  // We accumulate slave contributions to master DOFs and zero out slave
  // rows/cols.

  // Copy matrix to dense for easier manipulation (for small systems)
  // For large systems, a more efficient sparse manipulation would be needed
  Eigen::MatrixXcd A_dense = Eigen::MatrixXcd(asmbl.A);
  VecC b_new = asmbl.b;

  for (const auto &pair : pbc.pairs) {
    int m_edge = pair.master_edge;
    int s_edge = pair.slave_edge;
    double orient = static_cast<double>(pair.master_orient * pair.slave_orient);
    std::complex<double> phase = phi * orient;
    std::complex<double> phase_conj = phi_conj * orient;

    // Skip if either edge is in PEC BC (already constrained to zero)
    if (bc.dirichlet_edges.count(m_edge) || bc.dirichlet_edges.count(s_edge)) {
      continue;
    }

    // Accumulate slave row to master: A[m,:] += phase * A[s,:]
    A_dense.row(m_edge) += phase * A_dense.row(s_edge);

    // Accumulate slave column to master: A[:,m] += phase_conj * A[:,s]
    A_dense.col(m_edge) += phase_conj * A_dense.col(s_edge);

    // Accumulate RHS: b[m] += phase * b[s]
    b_new(m_edge) += phase * b_new(s_edge);
  }

  // Zero out slave rows and columns, set diagonal to 1
  for (const auto &pair : pbc.pairs) {
    int s_edge = pair.slave_edge;

    if (bc.dirichlet_edges.count(s_edge)) {
      continue;
    }

    // Zero row and column
    A_dense.row(s_edge).setZero();
    A_dense.col(s_edge).setZero();

    // Set diagonal to 1
    A_dense(s_edge, s_edge) = std::complex<double>(1.0, 0.0);

    // Set RHS to 0 (slave solution will be computed from master post-solve)
    b_new(s_edge) = std::complex<double>(0.0, 0.0);
  }

  // Convert back to sparse
  asmbl.A = A_dense.sparseView();
  asmbl.A.makeCompressed();
  asmbl.b = b_new;

  return asmbl;
}

Eigen::MatrixXcd
calculate_sparams_periodic(const Mesh &mesh, const MaxwellParams &p,
                           const BC &bc, const PeriodicBC &pbc,
                           const std::vector<WavePort> &ports) {
  const int num_ports = static_cast<int>(ports.size());
  Eigen::MatrixXcd S(num_ports, num_ports);

  // Build map of slave to master edges for solution recovery
  std::unordered_map<int, std::pair<int, std::complex<double>>> slave_to_master;
  for (const auto &pair : pbc.pairs) {
    double orient = static_cast<double>(pair.master_orient * pair.slave_orient);
    std::complex<double> phase = pbc.phase_shift * orient;
    slave_to_master[pair.slave_edge] = {pair.master_edge, phase};
  }

  for (int i = 0; i < num_ports; ++i) { // active port
    auto asmbl = assemble_maxwell_periodic(mesh, p, bc, pbc, ports, i);
    auto res = solve_linear(asmbl.A, asmbl.b, {});

    if (!check_solver_convergence(res,
                                  "calculate_sparams_periodic (active port " +
                                      std::to_string(i) + ")")) {
      for (int j = 0; j < num_ports; ++j) {
        S(j, i) =
            std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN());
      }
      continue;
    }

    // Recover slave DOF values from master DOFs
    VecC x_full = res.x;
    for (const auto &kv : slave_to_master) {
      int s_edge = kv.first;
      int m_edge = kv.second.first;
      std::complex<double> phase = kv.second.second;
      x_full(s_edge) = phase * x_full(m_edge);
    }

    fill_sparam_column(S, i, ports, bc, x_full);
  }
  return S;
}

Eigen::MatrixXcd
calculate_sparams_eigenmode(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, const std::vector<WavePort> &ports) {
  const int num_ports = static_cast<int>(ports.size());
  const int m = static_cast<int>(mesh.edges.size());
  Eigen::MatrixXcd S(num_ports, num_ports);

  if (num_ports == 0) {
    return S;
  }

  // Per-port propagation constants: β = sqrt(ε_r μ_r k0² − kc²) with the
  // material of the volume region adjacent to the port face (complex for
  // lossy media; principal sqrt gives Re>0, Im<0, i.e. a decaying
  // outgoing wave under the e^{+jωt} convention).
  const double k0 = p.omega / c0;
  std::vector<std::complex<double>> betas(num_ports);
  for (int i = 0; i < num_ports; ++i) {
    const double kc = ports[i].mode.kc;
    std::complex<double> eps_mu(1.0, 0.0);
    {
      // Find a tet with a face on the port surface to read its material.
      std::set<std::array<std::int64_t, 3>> port_faces;
      for (const auto &tri : mesh.tris) {
        if (tri.phys != ports[i].surface_tag)
          continue;
        std::array<std::int64_t, 3> f{tri.conn[0], tri.conn[1], tri.conn[2]};
        std::sort(f.begin(), f.end());
        port_faces.insert(f);
      }
      static const int tet_faces[4][3] = {
          {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
      bool found = false;
      for (const auto &tet : mesh.tets) {
        for (const auto &tf : tet_faces) {
          std::array<std::int64_t, 3> f{tet.conn[tf[0]], tet.conn[tf[1]],
                                        tet.conn[tf[2]]};
          std::sort(f.begin(), f.end());
          if (port_faces.count(f)) {
            eps_mu =
                p.get_eps_r(tet.phys, p.omega) * p.get_mu_r(tet.phys, p.omega);
            found = true;
            break;
          }
        }
        if (found)
          break;
      }
    }
    const std::complex<double> beta_sq = eps_mu * k0 * k0 - kc * kc;
    if (std::real(beta_sq) <= 0) {
      // Mode is evanescent at this frequency
      const double freq_hz = p.omega / (2.0 * M_PI);
      const double fc_hz =
          kc * c0 / (2.0 * M_PI * std::sqrt(std::abs(std::real(eps_mu))));
      std::cerr << "WARNING: Frequency " << freq_hz / 1e9
                << " GHz is below cutoff " << fc_hz / 1e9 << " GHz for port "
                << i
                << " — mode is evanescent, S-parameters may be meaningless."
                << std::endl;
      return S;
    }
    betas[i] = std::sqrt(beta_sq);
  }

  // Assembled port surface mass matrices M_s (global edge indexing).
  // The first-order modal ABC  n̂×(∇×E) + jβ E_t = 2jβ E_inc,t  contributes
  // jβ·M_s to the system matrix and 2jβ·M_s·e_inc to the RHS.
  std::vector<SpMatC> port_mass(num_ports);
  for (int i = 0; i < num_ports; ++i) {
    port_mass[i] = assemble_port_surface_mass(mesh, ports[i].surface_tag,
                                              bc.dirichlet_edges)
                       .cast<std::complex<double>>();
  }

  // Port eigenvectors normalized in the M_s inner product: e^H M_s e = 1,
  // so a unit-amplitude incident mode gives V_inc = 1 in the extraction.
  std::vector<Eigen::VectorXcd> port_vecs(num_ports);
  for (int i = 0; i < num_ports; ++i) {
    port_vecs[i] = Eigen::VectorXcd::Zero(m);
    for (size_t k = 0; k < ports[i].edges.size(); ++k) {
      int e = ports[i].edges[k];
      if (!bc.dirichlet_edges.count(e)) {
        port_vecs[i](e) = ports[i].weights(k);
      }
    }
    const double norm_sq =
        std::real(port_vecs[i].dot(port_mass[i] * port_vecs[i]));
    if (norm_sq > 1e-30) {
      port_vecs[i] /= std::sqrt(norm_sq);
    }
  }

  // Solve for each active port
  for (int active_port = 0; active_port < num_ports; ++active_port) {
    // Assemble base Maxwell system (no ports)
    std::vector<WavePort> no_ports;
    auto asmbl = assemble_maxwell(mesh, p, bc, no_ports, -1);

    // Add jβ·M_s on every port surface. port_abc_scale is 1.0 for the
    // assembled operator (kept as a diagnostic knob).
    for (int i = 0; i < num_ports; ++i) {
      const std::complex<double> coeff =
          std::complex<double>(0.0, p.port_abc_scale) * betas[i];
      for (int col = 0; col < port_mass[i].outerSize(); ++col) {
        for (SpMatC::InnerIterator it(port_mass[i], col); it; ++it) {
          asmbl.A.coeffRef(it.row(), it.col()) += coeff * it.value();
        }
      }
    }

    // Excitation: RHS = 2jβ·M_s·e_inc
    const std::complex<double> abc_active =
        std::complex<double>(0.0, p.port_abc_scale) * betas[active_port];
    asmbl.b =
        2.0 * abc_active * (port_mass[active_port] * port_vecs[active_port]);

    // Solve
    auto res = solve_linear(asmbl.A, asmbl.b, {});

    if (!check_solver_convergence(res,
                                  "calculate_sparams_eigenmode (active port " +
                                      std::to_string(active_port) + ")")) {
      for (int j = 0; j < num_ports; ++j) {
        S(j, active_port) =
            std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN());
      }
      continue;
    }

    // Extract modal amplitudes via M_s-weighted overlaps: with the
    // e^H M_s e = 1 normalization, Vj = e_j^H M_s x is the amplitude of
    // mode j in the solution and V_inc = 1 by construction.
    std::complex<double> V_inc(1.0, 0.0);

    for (int j = 0; j < num_ports; ++j) {
      std::complex<double> Vj = port_vecs[j].dot(port_mass[j] * res.x);

      if (j == active_port) {
        // Reflection: S_ii = (V_i - V_inc) / V_inc
        S(j, active_port) = (Vj - V_inc) / V_inc;
      } else {
        // Transmission: S_ji = V_j / V_inc
        S(j, active_port) = Vj / V_inc;
      }
    }
  }

  return S;
}

} // namespace edgefem
