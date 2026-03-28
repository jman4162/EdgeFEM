#include "edgefem/edge_basis.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <array>

namespace edgefem {
namespace {
const std::array<std::array<int, 2>, 6> edge_pairs = {
    std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2},
    std::array<int, 2>{0, 3}, std::array<int, 2>{1, 2},
    std::array<int, 2>{1, 3}, std::array<int, 2>{2, 3}};

void gradients_and_volume(const std::array<Eigen::Vector3d, 4> &v,
                          std::array<Eigen::Vector3d, 4> &g, double &V) {
  Eigen::Matrix3d B;
  B.col(0) = v[0] - v[3];
  B.col(1) = v[1] - v[3];
  B.col(2) = v[2] - v[3];
  Eigen::Matrix3d BinvT = B.inverse().transpose();
  g[0] = BinvT.col(0);
  g[1] = BinvT.col(1);
  g[2] = BinvT.col(2);
  g[3] = -g[0] - g[1] - g[2];
  V = std::abs(B.determinant()) / 6.0;
}

double lambda_int(int i, int j, double V) {
  return (i == j) ? V / 10.0 : V / 20.0;
}
} // namespace

Eigen::Matrix<double, 6, 3>
whitney_edge_curls(const std::array<Eigen::Vector3d, 4> &v) {
  std::array<Eigen::Vector3d, 4> g;
  double V;
  gradients_and_volume(v, g, V);
  Eigen::Matrix<double, 6, 3> curls;
  for (int i = 0; i < 6; ++i) {
    int a = edge_pairs[i][0];
    int b = edge_pairs[i][1];
    curls.row(i) = 2.0 * g[a].cross(g[b]);
  }
  return curls;
}

Eigen::Matrix<double, 6, 6>
whitney_curl_curl_matrix(const std::array<Eigen::Vector3d, 4> &v) {
  std::array<Eigen::Vector3d, 4> g;
  double V;
  gradients_and_volume(v, g, V);
  Eigen::Matrix<double, 6, 3> curls;
  for (int i = 0; i < 6; ++i) {
    int a = edge_pairs[i][0];
    int b = edge_pairs[i][1];
    curls.row(i) = 2.0 * g[a].cross(g[b]);
  }
  Eigen::Matrix<double, 6, 6> K;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      K(i, j) = V * curls.row(i).dot(curls.row(j));
  return K;
}

Eigen::Matrix<double, 6, 6>
whitney_mass_matrix(const std::array<Eigen::Vector3d, 4> &v) {
  std::array<Eigen::Vector3d, 4> g;
  double V;
  gradients_and_volume(v, g, V);
  Eigen::Matrix<double, 6, 6> M;
  for (int i = 0; i < 6; ++i) {
    int a = edge_pairs[i][0];
    int b = edge_pairs[i][1];
    for (int j = 0; j < 6; ++j) {
      int c = edge_pairs[j][0];
      int d = edge_pairs[j][1];
      double term = 0.0;
      term += g[b].dot(g[d]) * lambda_int(a, c, V);
      term -= g[b].dot(g[c]) * lambda_int(a, d, V);
      term -= g[a].dot(g[d]) * lambda_int(b, c, V);
      term += g[a].dot(g[c]) * lambda_int(b, d, V);
      M(i, j) = term;
    }
  }
  return M;
}

Eigen::Vector4d compute_barycentric(const std::array<Eigen::Vector3d, 4> &v,
                                    const Eigen::Vector3d &p) {
  // Build matrix T = [v0-v3, v1-v3, v2-v3]
  Eigen::Matrix3d T;
  T.col(0) = v[0] - v[3];
  T.col(1) = v[1] - v[3];
  T.col(2) = v[2] - v[3];

  // Solve T * [lambda0, lambda1, lambda2]^T = p - v3
  Eigen::Vector3d rhs = p - v[3];
  Eigen::Vector3d lambda012 = T.inverse() * rhs;

  Eigen::Vector4d lambda;
  lambda(0) = lambda012(0);
  lambda(1) = lambda012(1);
  lambda(2) = lambda012(2);
  lambda(3) = 1.0 - lambda(0) - lambda(1) - lambda(2);

  return lambda;
}

Eigen::Vector3d compute_grad_lambda(const std::array<Eigen::Vector3d, 4> &v,
                                    int i) {
  std::array<Eigen::Vector3d, 4> g;
  double V;
  gradients_and_volume(v, g, V);
  return g[i];
}

Eigen::Vector3cd
evaluate_edge_field(const std::array<Eigen::Vector3d, 4> &vertices,
                    const std::array<int, 6> &edge_orient,
                    const std::array<std::complex<double>, 6> &edge_dofs,
                    const Eigen::Vector3d &point) {

  Eigen::Vector4d lambda = compute_barycentric(vertices, point);

  Eigen::Vector3cd E = Eigen::Vector3cd::Zero();

  constexpr int edge_nodes[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                    {1, 2}, {1, 3}, {2, 3}};

  for (int e = 0; e < 6; ++e) {
    int n0 = edge_nodes[e][0];
    int n1 = edge_nodes[e][1];

    Eigen::Vector3d grad_l0 = compute_grad_lambda(vertices, n0);
    Eigen::Vector3d grad_l1 = compute_grad_lambda(vertices, n1);

    // Whitney basis: W_e = lambda_n0 * grad(lambda_n1) - lambda_n1 *
    // grad(lambda_n0)
    Eigen::Vector3d W_e = lambda(n0) * grad_l1 - lambda(n1) * grad_l0;

    std::complex<double> coeff =
        edge_dofs[e] * static_cast<double>(edge_orient[e]);

    E += coeff * W_e;
  }

  return E;
}

} // namespace edgefem
