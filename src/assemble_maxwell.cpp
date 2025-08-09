#include "vectorem/maxwell.hpp"
#include "vectorem/edge_basis.hpp"

#include <Eigen/SparseCore>
#include <vector>

namespace vectorem {

MaxwellAssembly assemble_maxwell(const Mesh &mesh,
                                 const MaxwellParams &p,
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
    auto Kloc = whitney_curl_curl_matrix(X);
    auto Mloc = whitney_mass_matrix(X);
    for (int i = 0; i < 6; ++i) {
      int gi = tet.edges[i];
      int si = tet.edge_orient[i];
      for (int j = 0; j < 6; ++j) {
        int gj = tet.edges[j];
        int sj = tet.edge_orient[j];
        std::complex<double> val =
            (Kloc(i, j) / p.mu_r) -
            (p.omega * p.omega * p.eps_r * Mloc(i, j));
        val *= static_cast<double>(si * sj);
        trips.emplace_back(gi, gj, val);
      }
    }
  }
  asmbl.A.setFromTriplets(trips.begin(), trips.end());

  for (int e : bc.dirichlet_edges) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A, e);
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
