#pragma once

#include <vector>

#include "edgefem/bc.hpp"
#include "edgefem/fem.hpp"

namespace edgefem {

/// Numbering of the edges that are not Dirichlet (PEC), for solving a
/// system with the Dirichlet rows and columns removed instead of kept as
/// identity rows.
///
/// With homogeneous Dirichlet conditions (E_t = 0 on PEC) the reduced
/// system A_ff x_f = b_f has the same solution on the free edges as the
/// full system with identity rows, and the full-length solution is x_f
/// scattered back with zeros on the Dirichlet edges. Free edges keep their
/// relative order.
class FreeDofs {
public:
  FreeDofs(int num_edges, const BC &bc);

  int num_full() const { return static_cast<int>(full_to_free_.size()); }
  int num_free() const { return static_cast<int>(free_to_full_.size()); }

  /// Position of edge e in the reduced system, or -1 for a Dirichlet edge.
  int free_index(int e) const { return full_to_free_[e]; }

  /// Rows and columns of the free edges of a num_full x num_full matrix.
  SpMatC restrict_matrix(const SpMatC &A) const;

  /// Entries of the free edges of a length num_full vector.
  VecC restrict_vector(const VecC &v) const;

  /// Length num_full vector: x_free on the free edges, zero elsewhere.
  VecC expand(const VecC &x_free) const;

private:
  std::vector<int> full_to_free_;
  std::vector<int> free_to_full_;
};

} // namespace edgefem
