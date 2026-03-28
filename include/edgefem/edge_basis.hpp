#pragma once

#include <Eigen/Core>
#include <array>

namespace edgefem {

Eigen::Matrix<double, 6, 3>
whitney_edge_curls(const std::array<Eigen::Vector3d, 4> &v);

Eigen::Matrix<double, 6, 6>
whitney_curl_curl_matrix(const std::array<Eigen::Vector3d, 4> &v);

Eigen::Matrix<double, 6, 6>
whitney_mass_matrix(const std::array<Eigen::Vector3d, 4> &v);

/// Compute barycentric coordinates of a point within a tetrahedron.
/// @param v The four vertices of the tetrahedron
/// @param p The point to compute coordinates for
/// @return Vector of 4 barycentric coordinates (lambda0, lambda1, lambda2,
/// lambda3)
Eigen::Vector4d compute_barycentric(const std::array<Eigen::Vector3d, 4> &v,
                                    const Eigen::Vector3d &p);

/// Compute gradient of barycentric coordinate lambda_i.
/// @param v The four vertices of the tetrahedron
/// @param i Index of the barycentric coordinate (0-3)
/// @return Gradient vector in 3D
Eigen::Vector3d compute_grad_lambda(const std::array<Eigen::Vector3d, 4> &v,
                                    int i);

/// Evaluate E-field at a point inside a tetrahedron from edge DOFs.
/// Uses Whitney (Nedelec) edge basis interpolation:
///   E(p) = sum_e coeff_e * W_e(p)
/// where W_e = lambda_n0 * grad(lambda_n1) - lambda_n1 * grad(lambda_n0).
///
/// @param vertices The four vertices of the tetrahedron
/// @param edge_orient Edge orientations (+1 or -1) for the 6 edges
/// @param edge_dofs Solution coefficients for the 6 edges
/// @param point The evaluation point (must be inside the tet)
/// @return Complex 3D electric field vector at the point
Eigen::Vector3cd
evaluate_edge_field(const std::array<Eigen::Vector3d, 4> &vertices,
                    const std::array<int, 6> &edge_orient,
                    const std::array<std::complex<double>, 6> &edge_dofs,
                    const Eigen::Vector3d &point);

} // namespace edgefem
