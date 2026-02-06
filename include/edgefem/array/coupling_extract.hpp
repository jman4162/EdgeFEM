#pragma once

#include <map>
#include <vector>

#include <Eigen/Core>

namespace edgefem {

/// Individual coupling coefficient between two ports.
struct CouplingCoefficient {
  int port_i;                       // First port index
  int port_j;                       // Second port index
  std::complex<double> S_ij;        // S-parameter value
  double magnitude_db;              // |S_ij| in dB
  double phase_deg;                 // arg(S_ij) in degrees
  double distance;                  // Distance between elements (meters)
};

/// Extract all coupling coefficients from S-matrix.
/// @param S N×N scattering matrix
/// @param positions Optional element positions for distance calculation
/// @return Vector of CouplingCoefficient for all pairs (i,j) where i != j
std::vector<CouplingCoefficient> extract_all_couplings(
    const Eigen::MatrixXcd &S,
    const std::vector<Eigen::Vector3d> &positions = {});

/// Compute coupling magnitude vs. element separation distance.
/// Useful for analyzing coupling decay with distance.
/// @param S N×N scattering matrix
/// @param positions Element positions (N x 3)
/// @return Map from distance (quantized to bins) to average coupling (dB)
std::map<int, double> coupling_vs_separation(
    const Eigen::MatrixXcd &S,
    const std::vector<Eigen::Vector3d> &positions);

/// Extract mutual coupling between nearest neighbors only.
/// @param S N×N scattering matrix
/// @param positions Element positions
/// @param max_neighbors Maximum number of neighbors to consider per element
/// @return Vector of CouplingCoefficient for nearest neighbors
std::vector<CouplingCoefficient> extract_nearest_neighbor_coupling(
    const Eigen::MatrixXcd &S,
    const std::vector<Eigen::Vector3d> &positions,
    int max_neighbors = 4);

/// Compute statistics on coupling levels.
struct CouplingStats {
  double max_coupling_db;        // Maximum coupling (excluding diagonal)
  double min_coupling_db;        // Minimum coupling
  double mean_coupling_db;       // Mean coupling
  double std_coupling_db;        // Standard deviation
  int max_coupling_i;            // Port index of max coupling (row)
  int max_coupling_j;            // Port index of max coupling (col)
};

/// Compute statistics on S-matrix coupling values.
/// @param S N×N scattering matrix
/// @return CouplingStats structure
CouplingStats compute_coupling_stats(const Eigen::MatrixXcd &S);

} // namespace edgefem
