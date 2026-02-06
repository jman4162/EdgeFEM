#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include "edgefem/array/active_impedance.hpp"
#include "edgefem/array/embedded_pattern.hpp"

namespace edgefem {

/// Complete array characterization package for export.
/// Compatible with Phased-Array-Antenna-Model import format.
struct ArrayPackage {
  std::string model_name;
  std::vector<Eigen::Vector3d> element_positions;
  std::vector<double> frequencies;
  std::vector<Eigen::MatrixXcd> S_matrices;  // One per frequency
  std::vector<std::vector<EmbeddedPattern>> embedded_patterns;  // [freq][port]
  std::vector<ActiveImpedanceResult> scan_data;
  double reference_impedance = 50.0;
};

/// Export complete array package to JSON format.
/// JSON format compatible with Phased-Array-Antenna-Model.
/// @param path Output file path
/// @param pkg Array package data
void export_array_package_json(const std::string &path, const ArrayPackage &pkg);

/// Export embedded element patterns to CSV format.
/// Columns: port_index, theta_deg, phi0_mag, phi0_phase, phi90_mag, phi90_phase
/// @param path Output file path
/// @param patterns Vector of embedded patterns
void export_patterns_csv(const std::string &path,
                         const std::vector<EmbeddedPattern> &patterns);

/// Export coupling matrix to CSV format.
/// Columns: port_i, port_j, S_real, S_imag, magnitude_db, phase_deg, distance
/// @param path Output file path
/// @param S Scattering matrix
/// @param positions Optional element positions
void export_coupling_csv(const std::string &path,
                         const Eigen::MatrixXcd &S,
                         const std::vector<Eigen::Vector3d> &positions = {});

/// Export scan angle sweep results to CSV.
/// Columns: theta_deg, phi_deg, element_idx, Gamma_real, Gamma_imag, Z_real, Z_imag, VSWR
/// @param path Output file path
/// @param scan_results Vector of scan results
void export_scan_csv(const std::string &path,
                     const std::vector<ActiveImpedanceResult> &scan_results);

} // namespace edgefem
