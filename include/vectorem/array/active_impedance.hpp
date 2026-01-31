#pragma once

#include <map>
#include <vector>

#include <Eigen/Core>

#include "vectorem/coupling.hpp"

namespace vectorem {

/// Result of active impedance calculation for one scan angle.
struct ActiveImpedanceResult {
  double theta;                    // Scan angle theta (radians)
  double phi;                      // Scan angle phi (radians)
  Eigen::VectorXcd Z_active;       // Active impedance per element
  Eigen::VectorXcd Gamma_active;   // Active reflection coefficient per element
  double scan_loss_db;             // Scan loss relative to broadside (dB)
};

/// Generate uniform phase excitation for scanning to (theta, phi).
/// @param num_elements Number of array elements
/// @param positions Element positions [N x 3]
/// @param theta Scan angle from broadside (radians)
/// @param phi Azimuth angle (radians)
/// @param k0 Free-space wavenumber (rad/m)
/// @return Complex excitation coefficients (N x 1)
Eigen::VectorXcd uniform_scan_excitation(int num_elements,
                                         const std::vector<Eigen::Vector3d> &positions,
                                         double theta, double phi, double k0);

/// Compute active impedance for given excitation.
/// @param coupling Coupling matrix (S, Z, Y parameters)
/// @param excitation Complex excitation coefficients
/// @return Active impedance result with Z, Gamma, and scan loss
ActiveImpedanceResult compute_active_impedance(const CouplingMatrix &coupling,
                                                const Eigen::VectorXcd &excitation);

/// Compute active impedance over a scan angle sweep.
/// @param coupling Coupling matrix at operating frequency
/// @param positions Element positions (meters)
/// @param theta_rad Vector of theta scan angles (radians)
/// @param phi_rad Vector of phi scan angles (radians)
/// @param k0 Free-space wavenumber (rad/m)
/// @return Vector of ActiveImpedanceResult, one per (theta, phi) pair
std::vector<ActiveImpedanceResult> active_impedance_scan(
    const CouplingMatrix &coupling,
    const std::vector<Eigen::Vector3d> &positions,
    const std::vector<double> &theta_rad,
    const std::vector<double> &phi_rad,
    double k0);

/// Compute active VSWR from active reflection coefficient.
/// @param Gamma_active Active reflection coefficient
/// @return VSWR (1.0 for perfect match, infinity for total reflection)
double active_vswr(std::complex<double> Gamma_active);

/// Find scan blindness angles where active VSWR exceeds threshold.
/// @param scan_results Vector of scan results
/// @param vswr_threshold VSWR threshold for blindness detection
/// @return Map of element index to vector of blindness angles
std::map<int, std::vector<double>> find_scan_blindness(
    const std::vector<ActiveImpedanceResult> &scan_results,
    double vswr_threshold = 10.0);

} // namespace vectorem
