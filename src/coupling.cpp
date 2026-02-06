#include "edgefem/coupling.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace edgefem {

CouplingMatrix compute_coupling(const Eigen::MatrixXcd &S, double freq,
                                double z0) {
  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }

  int n = static_cast<int>(S.rows());
  CouplingMatrix cm;
  cm.S = S;
  cm.frequency = freq;
  cm.z0 = z0;
  cm.num_ports = n;

  Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(n, n);

  // Z = z0 * (I + S) * (I - S)^-1
  Eigen::MatrixXcd I_minus_S = I - S;
  Eigen::MatrixXcd I_plus_S = I + S;

  // Use pseudoinverse for numerical stability
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXcd> cod(I_minus_S);
  if (cod.rank() < n) {
    // Matrix is singular, use regularization
    Eigen::MatrixXcd I_minus_S_reg = I_minus_S + 1e-12 * I;
    cm.Z = z0 * I_plus_S * I_minus_S_reg.inverse();
  } else {
    cm.Z = z0 * I_plus_S * cod.pseudoInverse();
  }

  // Y = (1/z0) * (I - S) * (I + S)^-1
  Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXcd> cod2(I_plus_S);
  if (cod2.rank() < n) {
    Eigen::MatrixXcd I_plus_S_reg = I_plus_S + 1e-12 * I;
    cm.Y = (1.0 / z0) * I_minus_S * I_plus_S_reg.inverse();
  } else {
    cm.Y = (1.0 / z0) * I_minus_S * cod2.pseudoInverse();
  }

  return cm;
}

Eigen::MatrixXd mutual_coupling_db(const CouplingMatrix &cm) {
  int n = cm.num_ports;
  Eigen::MatrixXd result(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double mag = std::abs(cm.S(i, j));
      // Avoid log(0) by clamping to a minimum value
      result(i, j) = 20.0 * std::log10(std::max(mag, 1e-20));
    }
  }

  return result;
}

Eigen::VectorXcd active_reflection(const Eigen::MatrixXcd &S,
                                   const Eigen::VectorXcd &excitation) {
  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }
  if (S.rows() != excitation.size()) {
    throw std::runtime_error(
        "Excitation vector size must match S-matrix dimension");
  }

  int n = static_cast<int>(S.rows());
  Eigen::VectorXcd gamma(n);

  // Gamma_active[i] = sum_j(S[i,j] * a[j]) / a[i]
  // where a[j] is the excitation coefficient for port j
  for (int i = 0; i < n; ++i) {
    std::complex<double> numerator = 0.0;
    for (int j = 0; j < n; ++j) {
      numerator += S(i, j) * excitation(j);
    }
    if (std::abs(excitation(i)) < 1e-20) {
      // Port i is not excited, set reflection to 0
      gamma(i) = 0.0;
    } else {
      gamma(i) = numerator / excitation(i);
    }
  }

  return gamma;
}

bool is_passive(const Eigen::MatrixXcd &S, double tolerance) {
  if (S.rows() != S.cols()) {
    return false;
  }

  // A network is passive if all singular values of S are <= 1
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(S);
  Eigen::VectorXd singular_values = svd.singularValues();

  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values(i) > 1.0 + tolerance) {
      return false;
    }
  }

  return true;
}

bool is_lossless(const Eigen::MatrixXcd &S, double tolerance) {
  if (S.rows() != S.cols()) {
    return false;
  }

  // A network is lossless if S^H * S = I
  int n = static_cast<int>(S.rows());
  Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(n, n);
  Eigen::MatrixXcd product = S.adjoint() * S;

  return product.isApprox(I, tolerance);
}

bool is_reciprocal(const Eigen::MatrixXcd &S, double tolerance) {
  if (S.rows() != S.cols()) {
    return false;
  }

  // A network is reciprocal if S = S^T
  return S.isApprox(S.transpose(), tolerance);
}

} // namespace edgefem
