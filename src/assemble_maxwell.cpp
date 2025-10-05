#include "vectorem/edge_basis.hpp"
#include "vectorem/maxwell.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "vectorem/solver.hpp"

namespace vectorem {

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
    Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d max = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
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
      const double integral = sigma_max * thickness / (kv.second.grading_order + 1.0);
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
        Eigen::Vector3d centroid =
            (X[0] + X[1] + X[2] + X[3]) / 4.0;
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
          if (p.enforce_pml_heuristics && sigma_max > 0.0 && sigma[axis] < sigma_max * kMinProfile) {
            sigma[axis] = sigma_max * kMinProfile;
          }
          stretch_tensor[axis] = std::complex<double>(1.0, sigma[axis] / p.omega);
        }
        std::complex<double> stretch_sum = stretch_tensor[0] + stretch_tensor[1] + stretch_tensor[2];
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
    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int si = tet.edge_orient[i];
      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int sj = tet.edge_orient[j];
        std::complex<double> val =
            (Kloc(i, j) / p.mu_r) - (p.omega * p.omega * p.eps_r * Mloc(i, j));
        val *= static_cast<double>(si * sj);
        trips.emplace_back(gi, gj, val);
      }
    }
  }
  asmbl.A.setFromTriplets(trips.begin(), trips.end());

  asmbl.A.makeCompressed();

  // Add port admittances and sources
  for (size_t i = 0; i < ports.size(); ++i) {
    const auto &port = ports[i];
    if (port.weights.size() != static_cast<int>(port.edges.size()))
      continue;
    if (port.mode.Z0 == std::complex<double>(0.0))
      continue;

    for (size_t a = 0; a < port.edges.size(); ++a) {
      int ga = port.edges[a];
      const std::complex<double> wa = port.weights[a];
      for (size_t b = 0; b < port.edges.size(); ++b) {
        int gb = port.edges[b];
        const std::complex<double> wb = port.weights[b];
        asmbl.A.coeffRef(ga, gb) += wa * std::conj(wb) / port.mode.Z0;
      }
    }

    if (static_cast<int>(i) == active_port_idx) {
      const std::complex<double> scale = 2.0 / std::sqrt(port.mode.Z0);
      for (size_t a = 0; a < port.edges.size(); ++a) {
        int ga = port.edges[a];
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

  for (int e : bc.dirichlet_edges) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A,
                                                                     e);
         it; ++it) {
      it.valueRef() = (it.row() == e && it.col() == e)
                          ? std::complex<double>(1.0, 0.0)
                          : std::complex<double>(0.0, 0.0);
    }
    asmbl.b(e) = 0.0;
  }
  return asmbl;
}

Eigen::MatrixXcd
calculate_sparams(const Mesh &mesh, const MaxwellParams &p, const BC &bc,
                  const std::vector<WavePort> &ports) {
  const int num_ports = ports.size();
  Eigen::MatrixXcd S(num_ports, num_ports);

  for (int i = 0; i < num_ports; ++i) { // active port
    auto asmbl = assemble_maxwell(mesh, p, bc, ports, i);
    auto res = solve_linear(asmbl.A, asmbl.b, {});

    const std::complex<double> V_inc_i = std::sqrt(ports[i].mode.Z0);

    for (int j = 0; j < num_ports; ++j) { // passive port
      std::complex<double> Vj = 0.0;
      for (size_t k = 0; k < ports[j].edges.size(); ++k) {
        int edge_idx = ports[j].edges[k];
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

} // namespace vectorem
