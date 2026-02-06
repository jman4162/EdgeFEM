#pragma once

#include <Eigen/Core>
#include <complex>

namespace edgefem {

/// Coupling matrix representation in multiple network parameter formats.
struct CouplingMatrix {
  Eigen::MatrixXcd Z; // Impedance matrix
  Eigen::MatrixXcd Y; // Admittance matrix
  Eigen::MatrixXcd S; // Scattering matrix
  double frequency;   // Operating frequency (Hz)
  double z0;          // Reference impedance (ohms)
  int num_ports;      // Number of ports
};

/// Compute impedance/admittance matrices from S-parameters.
/// @param S N×N scattering matrix
/// @param freq Operating frequency in Hz
/// @param z0 Reference impedance (default 50 ohms)
/// @return CouplingMatrix with Z, Y, S matrices
CouplingMatrix compute_coupling(const Eigen::MatrixXcd &S, double freq,
                                double z0 = 50.0);

/// Convert S-parameter matrix to mutual coupling in dB.
/// Returns a matrix where element (i,j) = 20*log10(|S(i,j)|).
/// @param cm CouplingMatrix containing S-parameters
/// @return N×N matrix of coupling values in dB
Eigen::MatrixXd mutual_coupling_db(const CouplingMatrix &cm);

/// Compute active reflection coefficient for each port given an excitation.
/// Active reflection accounts for all ports being excited simultaneously.
/// @param S N×N scattering matrix
/// @param excitation N×1 vector of port excitation coefficients
/// @return N×1 vector of active reflection coefficients
Eigen::VectorXcd active_reflection(const Eigen::MatrixXcd &S,
                                   const Eigen::VectorXcd &excitation);

/// Check if a network is passive.
/// A network is passive if |S| eigenvalues are all <= 1.
/// @param S N×N scattering matrix
/// @param tolerance Numerical tolerance for passivity check
/// @return true if network is passive, false otherwise
bool is_passive(const Eigen::MatrixXcd &S, double tolerance = 1e-6);

/// Check if a network is lossless.
/// A network is lossless if S^H * S = I.
/// @param S N×N scattering matrix
/// @param tolerance Numerical tolerance for lossless check
/// @return true if network is lossless, false otherwise
bool is_lossless(const Eigen::MatrixXcd &S, double tolerance = 1e-6);

/// Check if a network is reciprocal.
/// A network is reciprocal if S = S^T.
/// @param S N×N scattering matrix
/// @param tolerance Numerical tolerance for reciprocity check
/// @return true if network is reciprocal, false otherwise
bool is_reciprocal(const Eigen::MatrixXcd &S, double tolerance = 1e-6);

} // namespace edgefem
