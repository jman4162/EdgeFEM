#include "edgefem/edge_basis.hpp"
#include "edgefem/maxwell.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "edgefem/solver.hpp"

namespace edgefem {

namespace {
constexpr double c0 = 299792458.0;        // speed of light in vacuum (m/s)
constexpr double mu0 = 4.0 * M_PI * 1e-7; // permeability of free space (H/m)
constexpr double eps0 =
    1.0 / (mu0 * c0 * c0);        // permittivity of free space (F/m)
constexpr double eta0 = mu0 * c0; // free-space impedance ≈ 377 ohms
} // namespace

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

  // Add port admittances and sources
  // NOTE: Skip edges that are constrained by Dirichlet BC (PEC) to avoid
  // inconsistent system. PEC edges have zero field by definition.
  for (size_t i = 0; i < ports.size(); ++i) {
    const auto &port = ports[i];
    if (port.weights.size() != static_cast<int>(port.edges.size()))
      continue;
    if (port.mode.Z0 == std::complex<double>(0.0))
      continue;

    for (size_t a = 0; a < port.edges.size(); ++a) {
      int ga = port.edges[a];
      if (bc.dirichlet_edges.count(ga))
        continue; // Skip PEC edges
      const std::complex<double> wa = port.weights[a];
      for (size_t b = 0; b < port.edges.size(); ++b) {
        int gb = port.edges[b];
        if (bc.dirichlet_edges.count(gb))
          continue; // Skip PEC edges
        const std::complex<double> wb = port.weights[b];
        asmbl.A.coeffRef(ga, gb) += wa * std::conj(wb) / port.mode.Z0;
      }
    }

    if (static_cast<int>(i) == active_port_idx) {
      const std::complex<double> scale = 2.0 / std::sqrt(port.mode.Z0);
      for (size_t a = 0; a < port.edges.size(); ++a) {
        int ga = port.edges[a];
        if (bc.dirichlet_edges.count(ga))
          continue; // Skip PEC edges
        asmbl.b(ga) += scale * port.weights[a];
      }
    }
  }

  if (p.use_abc) {
    // Simple first-order absorbing boundary: add imaginary damping to
    // diagonal entries of boundary edges not marked as PEC.
    std::unordered_set<int> boundary_edges;
    for (const auto &tri : mesh.tris) {
      for (int e = 0; e < 3; ++e) {
        int ge = tri.edges[e];
        if (!bc.dirichlet_edges.count(ge))
          boundary_edges.insert(ge);
      }
    }
    for (int e : boundary_edges) {
      asmbl.A.coeffRef(e, e) += std::complex<double>(0.0, p.omega);
    }
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

Eigen::MatrixXcd calculate_sparams(const Mesh &mesh, const MaxwellParams &p,
                                   const BC &bc,
                                   const std::vector<WavePort> &ports) {
  const int num_ports = static_cast<int>(ports.size());
  Eigen::MatrixXcd S(num_ports, num_ports);

  for (int i = 0; i < num_ports; ++i) { // active port
    auto asmbl = assemble_maxwell(mesh, p, bc, ports, i);
    auto res = solve_linear(asmbl.A, asmbl.b, {});

    const std::complex<double> V_inc_i = std::sqrt(ports[i].mode.Z0);

    for (int j = 0; j < num_ports; ++j) { // passive port
      std::complex<double> Vj = 0.0;
      for (size_t k = 0; k < ports[j].edges.size(); ++k) {
        int edge_idx = ports[j].edges[k];
        // Skip PEC edges (they have zero solution by definition)
        if (bc.dirichlet_edges.count(edge_idx))
          continue;
        Vj += std::conj(ports[j].weights[k]) * res.x(edge_idx);
      }
      if (i == j) { // Reflection coefficient S_ii
        const std::complex<double> V_ref_i = Vj - V_inc_i;
        S(j, i) = V_ref_i / V_inc_i;
      } else { // Transmission coefficient S_ji
        S(j, i) = Vj / V_inc_i;
      }
    }
  }
  return S;
}

void normalize_port_weights(const Mesh &mesh, const MaxwellParams &p,
                            const BC &bc, std::vector<WavePort> &ports) {
  if (ports.empty())
    return;

  // Assemble FEM matrix without port loading
  std::vector<WavePort> no_ports;
  auto asmbl = assemble_maxwell(mesh, p, bc, no_ports, -1);

  for (auto &port : ports) {
    if (port.weights.size() != static_cast<int>(port.edges.size()))
      continue;
    if (port.mode.Z0 == std::complex<double>(0.0))
      continue;

    const int n_edges = static_cast<int>(port.edges.size());
    const double Z0_real = std::real(port.mode.Z0);

    // Compute w^H A^-1 w exactly via forward solve
    // Construct RHS = w (only at port edges)
    const int m = asmbl.A.rows();
    VecC rhs = VecC::Zero(m);
    for (int k = 0; k < n_edges; ++k) {
      int edge_idx = port.edges[k];
      if (!bc.dirichlet_edges.count(edge_idx)) {
        rhs(edge_idx) = port.weights[k];
      }
    }

    // Solve A * y = w to get y = A^-1 * w
    auto sol = solve_linear(asmbl.A, rhs, {});

    // Compute w^H * y = w^H * A^-1 * w
    std::complex<double> wAinvw = 0.0;
    int free_count = 0;
    for (int k = 0; k < n_edges; ++k) {
      int edge_idx = port.edges[k];
      if (!bc.dirichlet_edges.count(edge_idx)) {
        wAinvw += std::conj(port.weights[k]) * sol.x(edge_idx);
        free_count++;
      }
    }

    double wAinvw_mag = std::abs(wAinvw);

    // Also compute diagonal approximation for comparison
    double wAinvw_diag_est = 0.0;
    for (int k = 0; k < n_edges; ++k) {
      int edge_idx = port.edges[k];
      if (!bc.dirichlet_edges.count(edge_idx)) {
        std::complex<double> w_k = port.weights[k];
        double w_sq = std::real(std::conj(w_k) * w_k);
        double A_diag = std::abs(asmbl.A.coeff(edge_idx, edge_idx));
        if (A_diag > 1e-15) {
          wAinvw_diag_est += w_sq / A_diag;
        }
      }
    }

    // Target: w^H A^-1 w = Z0 for proper S-parameter extraction
    // This comes from the Woodbury identity analysis:
    // With port loading Y = ww^H/Z0, using b = 2w/sqrt(Z0), the port voltage
    // is: V = w^H (A+Y)^-1 b = 2*Z0*w^H A^-1 w / (sqrt(Z0)*(Z0 + w^H A^-1 w))
    // For matched port (V = V_inc = sqrt(Z0)), we need w^H A^-1 w = Z0
    double target = Z0_real;
    double ratio = wAinvw_mag / target;

    std::cerr << "  Port normalization: free_edges=" << free_count
              << ", wAinvw_exact=" << wAinvw_mag
              << ", wAinvw_diag=" << wAinvw_diag_est
              << ", target=Z0/2=" << target << ", ratio=" << ratio << std::endl;

    if (wAinvw_mag > 1e-15 && target > 1e-15) {
      // Scale factor: alpha^2 * wAinvw = target => alpha = sqrt(target/wAinvw)
      double alpha = std::sqrt(target / wAinvw_mag);
      std::cerr << "  Applying scale factor alpha=" << alpha << std::endl;
      port.weights *= alpha;
    }
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

  const int m = static_cast<int>(mesh.edges.size());
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

    // Recover slave DOF values from master DOFs
    VecC x_full = res.x;
    for (const auto &kv : slave_to_master) {
      int s_edge = kv.first;
      int m_edge = kv.second.first;
      std::complex<double> phase = kv.second.second;
      x_full(s_edge) = phase * x_full(m_edge);
    }

    const std::complex<double> V_inc_i = std::sqrt(ports[i].mode.Z0);

    for (int j = 0; j < num_ports; ++j) { // passive port
      std::complex<double> Vj = 0.0;
      for (size_t k = 0; k < ports[j].edges.size(); ++k) {
        int edge_idx = ports[j].edges[k];
        // Skip PEC edges (they have zero solution by definition)
        if (bc.dirichlet_edges.count(edge_idx))
          continue;
        Vj += std::conj(ports[j].weights[k]) * x_full(edge_idx);
      }
      if (i == j) {
        const std::complex<double> V_ref_i = Vj - V_inc_i;
        S(j, i) = V_ref_i / V_inc_i;
      } else {
        S(j, i) = Vj / V_inc_i;
      }
    }
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

  // Compute propagation constant from first port's mode parameters
  const double k0 = p.omega / c0;
  const double kc = ports[0].mode.kc;
  const double beta_sq = k0 * k0 - kc * kc;

  if (beta_sq <= 0) {
    // Mode is evanescent at this frequency
    return S;
  }
  const double beta = std::sqrt(beta_sq);

  // Create normalized port eigenvectors from port weights
  std::vector<Eigen::VectorXcd> port_vecs(num_ports);

  for (int i = 0; i < num_ports; ++i) {
    port_vecs[i] = Eigen::VectorXcd::Zero(m);

    // Copy weights to full-size vector
    for (size_t k = 0; k < ports[i].edges.size(); ++k) {
      int e = ports[i].edges[k];
      if (!bc.dirichlet_edges.count(e)) {
        port_vecs[i](e) = ports[i].weights(k);
      }
    }

    // Normalize
    double norm = port_vecs[i].norm();
    if (norm > 1e-15) {
      port_vecs[i] /= norm;
    }
  }

  // ABC coefficient with optimal scale factor (0.5 gives best results)
  const std::complex<double> abc_coeff(0.0, beta * p.port_abc_scale);

  // Solve for each active port
  for (int active_port = 0; active_port < num_ports; ++active_port) {
    // Assemble base Maxwell system (no ports)
    std::vector<WavePort> no_ports;
    auto asmbl = assemble_maxwell(mesh, p, bc, no_ports, -1);

    // Add ABC on all port edges
    for (int i = 0; i < num_ports; ++i) {
      for (int e : ports[i].edges) {
        if (!bc.dirichlet_edges.count(e)) {
          asmbl.A.coeffRef(e, e) += abc_coeff;
        }
      }
    }

    // Excitation: RHS = 2*abc_coeff * v_active
    asmbl.b = 2.0 * abc_coeff * port_vecs[active_port];

    // Solve
    auto res = solve_linear(asmbl.A, asmbl.b, {});

    // Extract S-parameters from overlap integrals
    std::complex<double> V_inc(
        1.0, 0.0); // Unit incident amplitude (normalized eigenvector)

    for (int j = 0; j < num_ports; ++j) {
      std::complex<double> Vj = port_vecs[j].dot(res.x);

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
