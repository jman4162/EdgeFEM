#include "vectorem/edge_basis.hpp"
#include "vectorem/maxwell.hpp"

#include <Eigen/SparseCore>
#include <unordered_set>
#include <vector>

namespace vectorem {

MaxwellAssembly assemble_maxwell(const Mesh &mesh, const MaxwellParams &p,
                                 const BC &bc) {
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

} // namespace vectorem
