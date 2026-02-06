#include "edgefem/fem.hpp"

#include <Eigen/SparseCore>

namespace edgefem {

ScalarAssembly assemble_scalar_helmholtz(const Mesh &mesh,
                                         const ScalarHelmholtzParams &p,
                                         const BC &bc) {
  const int n = static_cast<int>(mesh.nodes.size());
  ScalarAssembly asmbl;
  asmbl.A.resize(n, n);
  asmbl.b = VecC::Zero(n);

  // Placeholder: identity matrix
  std::vector<Eigen::Triplet<std::complex<double>>> trips;
  trips.reserve(n);
  for (int i = 0; i < n; ++i) {
    trips.emplace_back(i, i, std::complex<double>(1.0, 0.0));
  }
  asmbl.A.setFromTriplets(trips.begin(), trips.end());

  // Apply Dirichlet conditions by zeroing rows and cols, setting diagonal to 1
  for (int idx : bc.dirichlet_nodes) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(asmbl.A,
                                                                     idx);
         it; ++it) {
      it.valueRef() = (it.row() == idx && it.col() == idx)
                          ? std::complex<double>(1.0, 0.0)
                          : std::complex<double>(0.0, 0.0);
    }
    asmbl.b(idx) = 0.0;
  }
  return asmbl;
}

} // namespace edgefem
