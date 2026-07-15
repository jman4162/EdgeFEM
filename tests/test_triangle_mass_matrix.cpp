// Unit tests for triangle_whitney_mass_matrix.
//
// Checks the closed-form 3x3 edge-element mass matrix on a triangle against
// numerical quadrature of ∫ N_e · N_f dS, plus symmetry and positive
// definiteness. Whitney basis N_e = λ_a ∇λ_b − λ_b ∇λ_a with local edges
// (0,1), (1,2), (2,0) — the ordering build_edges() uses for Tri3.

#include "edgefem/edge_basis.hpp"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace edgefem;

namespace {

// Evaluate all three Whitney edge basis functions at barycentric (l0,l1,l2).
std::array<Eigen::Vector3d, 3>
whitney_at(const std::array<Eigen::Vector3d, 3> &v, double l0, double l1,
           double l2) {
  Eigen::Vector3d normal = (v[1] - v[0]).cross(v[2] - v[0]);
  double area2 = normal.norm();
  Eigen::Vector3d n_hat = normal / area2;
  std::array<Eigen::Vector3d, 3> g;
  g[0] = n_hat.cross(v[2] - v[1]) / area2;
  g[1] = n_hat.cross(v[0] - v[2]) / area2;
  g[2] = n_hat.cross(v[1] - v[0]) / area2;
  double lam[3] = {l0, l1, l2};
  static const int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  std::array<Eigen::Vector3d, 3> N;
  for (int e = 0; e < 3; ++e) {
    int a = pairs[e][0], b = pairs[e][1];
    N[e] = lam[a] * g[b] - lam[b] * g[a];
  }
  return N;
}

// Numerical quadrature of M[e][f] = ∫ N_e·N_f dS using a 7-point degree-5
// rule (exact for the degree-2 integrand).
Eigen::Matrix3d quadrature_mass(const std::array<Eigen::Vector3d, 3> &v) {
  const double w0 = 9.0 / 40.0;
  const double a1 = 0.059715871789770, b1 = 0.470142064105115,
               w1 = 0.132394152788506;
  const double a2 = 0.797426985353087, b2 = 0.101286507323456,
               w2 = 0.125939180544827;
  const double pts[7][4] = {{1.0 / 3, 1.0 / 3, 1.0 / 3, w0},
                            {a1, b1, b1, w1},
                            {b1, a1, b1, w1},
                            {b1, b1, a1, w1},
                            {a2, b2, b2, w2},
                            {b2, a2, b2, w2},
                            {b2, b2, a2, w2}};
  double area = 0.5 * (v[1] - v[0]).cross(v[2] - v[0]).norm();
  Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
  for (const auto &p : pts) {
    auto N = whitney_at(v, p[0], p[1], p[2]);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        M(i, j) += p[3] * area * N[i].dot(N[j]);
  }
  return M;
}

void check_triangle(const std::array<Eigen::Vector3d, 3> &v,
                    const char *label) {
  Eigen::Matrix3d M = triangle_whitney_mass_matrix(v);
  Eigen::Matrix3d Mq = quadrature_mass(v);

  double err = (M - Mq).norm() / Mq.norm();
  std::cout << label << ": closed-form vs quadrature rel. error = " << err
            << std::endl;
  assert(err < 1e-12);

  assert((M - M.transpose()).norm() < 1e-14 * M.norm());

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M);
  double min_ev = es.eigenvalues().minCoeff();
  std::cout << label << ": min eigenvalue = " << min_ev << std::endl;
  assert(min_ev > 0.0); // 3 tangential edge functions are independent
}

} // namespace

int main() {
  // Reference right triangle in the xy-plane
  check_triangle({Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0),
                  Eigen::Vector3d(0, 1, 0)},
                 "unit right triangle");

  // Scaled, shifted, non-right triangle
  check_triangle({Eigen::Vector3d(0.2, -0.1, 0.0),
                  Eigen::Vector3d(1.7, 0.3, 0.0),
                  Eigen::Vector3d(0.4, 2.1, 0.0)},
                 "general planar triangle");

  // Triangle embedded with arbitrary 3D orientation
  check_triangle({Eigen::Vector3d(0.1, 0.2, 0.3),
                  Eigen::Vector3d(1.1, 0.4, 0.9),
                  Eigen::Vector3d(0.3, 1.5, 0.7)},
                 "3D-embedded triangle");

  std::cout << "All triangle mass matrix tests passed." << std::endl;
  return 0;
}
