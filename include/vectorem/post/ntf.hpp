#pragma once

#include <complex>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace vectorem {

/// Far-field point from near-to-far transformation (Stratton-Chu).
/// Contains E_theta and E_phi far-field components.
struct NTFPoint2D {
  double theta_deg;         ///< Polar angle in degrees
  std::complex<double> e_theta; ///< \f$E_\theta\f$ far-field component
  std::complex<double> e_phi;   ///< \f$E_\phi\f$ far-field component
};

/// 3D far-field pattern over spherical (theta, phi) grid.
/// Matrices are Ntheta x Nphi where theta varies along rows, phi along columns.
struct FFPattern3D {
  Eigen::MatrixXd theta_grid;    ///< Theta angles in radians
  Eigen::MatrixXd phi_grid;      ///< Phi angles in radians
  Eigen::MatrixXcd E_theta;      ///< E_theta far-field component
  Eigen::MatrixXcd E_phi;        ///< E_phi far-field component

  /// Get total electric field magnitude at each point
  Eigen::MatrixXd total_magnitude() const;

  /// Get total power pattern (|E_theta|^2 + |E_phi|^2)
  Eigen::MatrixXd power_pattern() const;

  /// Get pattern in dB relative to maximum (normalized)
  Eigen::MatrixXd pattern_dB() const;
};

/**
 * Compute far-field pattern using a discretized Stratton--Chu integral over a
 * Huygens surface. The surface is represented by sample points with outward
 * normals, tangential electric and magnetic fields, and associated patch
 * areas.
 *
 * The returned pattern is sampled for a fixed azimuthal angle (phi) and a list
 * of polar angles (theta).
 */
std::vector<NTFPoint2D> stratton_chu_2d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    double phi_rad, double k0);

/**
 * Compute full 3D far-field pattern using Stratton-Chu integral over a
 * Huygens surface. Returns pattern sampled over spherical (theta, phi) grid.
 *
 * @param r      Surface sample point positions
 * @param n      Outward surface normals at each point
 * @param E      Tangential electric field at each point
 * @param H      Tangential magnetic field at each point
 * @param area   Surface patch area at each point
 * @param theta_rad  Vector of theta angles to sample (0 to pi)
 * @param phi_rad    Vector of phi angles to sample (0 to 2*pi)
 * @param k0     Free-space wavenumber (2*pi/lambda)
 * @return FFPattern3D with Ntheta x Nphi grids
 */
FFPattern3D stratton_chu_3d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    const std::vector<double> &phi_rad,
    double k0);

/**
 * Compute directivity from 3D far-field pattern.
 * Directivity D = 4*pi * U_max / P_rad
 * where U is radiation intensity and P_rad is total radiated power.
 *
 * @param pattern  3D far-field pattern
 * @return Directivity (dimensionless, linear scale)
 */
double compute_directivity(const FFPattern3D &pattern);

/**
 * Compute maximum gain from 3D pattern with radiation efficiency.
 * G_max = efficiency * D_max
 *
 * @param pattern    3D far-field pattern
 * @param efficiency Radiation efficiency (0 to 1, default 1.0)
 * @return Maximum gain (dimensionless, linear scale)
 */
double compute_max_gain(const FFPattern3D &pattern, double efficiency = 1.0);

/**
 * Compute half-power beamwidth (HPBW) in E-plane and H-plane.
 *
 * @param pattern 3D far-field pattern
 * @return Pair of (E-plane HPBW, H-plane HPBW) in degrees
 */
std::pair<double, double> compute_hpbw(const FFPattern3D &pattern);

/**
 * Export 3D pattern to CSV file with columns:
 * theta_deg, phi_deg, E_theta_mag, E_theta_phase, E_phi_mag, E_phi_phase, total_dB
 */
void write_pattern_3d_csv(const std::string &path, const FFPattern3D &pattern);

/**
 * Export 3D pattern to VTK file for ParaView visualization.
 * Creates a spherical surface colored by pattern magnitude.
 */
void write_pattern_3d_vtk(const std::string &path, const FFPattern3D &pattern,
                          double radius = 1.0);

/** Write a 2D pattern to CSV for visualization. */
void write_pattern_csv(const std::string &path, double phi_deg,
                       const std::vector<NTFPoint2D> &pat);

} // namespace vectorem

