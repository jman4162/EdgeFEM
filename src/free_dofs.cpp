#include "edgefem/free_dofs.hpp"

#include <stdexcept>
#include <string>

namespace edgefem {

FreeDofs::FreeDofs(int num_edges, const BC &bc)
    : full_to_free_(static_cast<std::size_t>(num_edges), -1) {
  free_to_full_.reserve(static_cast<std::size_t>(num_edges));
  for (int e = 0; e < num_edges; ++e) {
    if (bc.dirichlet_edges.count(e))
      continue;
    full_to_free_[e] = static_cast<int>(free_to_full_.size());
    free_to_full_.push_back(e);
  }
}

SpMatC FreeDofs::restrict_matrix(const SpMatC &A) const {
  if (A.rows() != num_full() || A.cols() != num_full())
    throw std::invalid_argument(
        "FreeDofs::restrict_matrix: matrix is " + std::to_string(A.rows()) +
        " x " + std::to_string(A.cols()) + ", expected " +
        std::to_string(num_full()) + " x " + std::to_string(num_full()));
  const int n = num_free();
  SpMatC out(n, n);
  // Column-major: copy each free column, dropping Dirichlet rows. Row order
  // within a column is preserved, so the result is already sorted.
  Eigen::VectorXi nnz_per_col(n);
  for (int jf = 0; jf < n; ++jf) {
    int count = 0;
    for (SpMatC::InnerIterator it(A, free_to_full_[jf]); it; ++it)
      if (full_to_free_[it.row()] >= 0)
        ++count;
    nnz_per_col(jf) = count;
  }
  out.reserve(nnz_per_col);
  for (int jf = 0; jf < n; ++jf)
    for (SpMatC::InnerIterator it(A, free_to_full_[jf]); it; ++it) {
      const int i = full_to_free_[it.row()];
      if (i >= 0)
        out.insertBackUncompressed(i, jf) = it.value();
    }
  out.makeCompressed();
  return out;
}

VecC FreeDofs::restrict_vector(const VecC &v) const {
  if (v.size() != num_full())
    throw std::invalid_argument(
        "FreeDofs::restrict_vector: vector has " + std::to_string(v.size()) +
        " entries, expected " + std::to_string(num_full()));
  VecC out(num_free());
  for (int i = 0; i < num_free(); ++i)
    out(i) = v(free_to_full_[i]);
  return out;
}

VecC FreeDofs::expand(const VecC &x_free) const {
  if (x_free.size() != num_free())
    throw std::invalid_argument(
        "FreeDofs::expand: vector has " + std::to_string(x_free.size()) +
        " entries, expected " + std::to_string(num_free()));
  VecC out = VecC::Zero(num_full());
  for (int i = 0; i < num_free(); ++i)
    out(free_to_full_[i]) = x_free(i);
  return out;
}

} // namespace edgefem
