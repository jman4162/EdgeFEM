#include "vectorem/edge_basis.hpp"
#include "vectorem/maxwell.hpp"

#include <Eigen/SparseCore>
#include <complex>
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

  for (const auto &tet : mesh.tets) {
    std::array<Eigen::Vector3d, 4> X;
    for (int i = 0; i < 4; ++i) {
      int idx = mesh.nodeIndex.at(tet.conn[i]);
      X[i] = mesh.nodes[idx].xyz;
    }
    // Determine if this element lies in a PML region
    std::complex<double> stretch(1.0, 0.0);
    if (p.omega != 0.0 && p.pml_regions.count(tet.phys)) {
      stretch = std::complex<double>(1.0, p.pml_sigma / p.omega);
    }
    Eigen::Matrix<std::complex<double>, 6, 6> Kloc =
        whitney_curl_curl_matrix(X).cast<std::complex<double>>() / stretch;
    Eigen::Matrix<std::complex<double>, 6, 6> Mloc =
        whitney_mass_matrix(X).cast<std::complex<double>>() * stretch;
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
