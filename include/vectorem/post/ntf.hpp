#pragma once

#include <complex>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace vectorem {

struct FFPoint2D {
  double theta_deg;         ///< Polar angle in degrees
  std::complex<double> e_theta; ///< \f$E_\theta\f$ far-field component
  std::complex<double> e_phi;   ///< \f$E_\phi\f$ far-field component
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
std::vector<FFPoint2D> stratton_chu_2d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    double phi_rad, double k0);

/** Write a 2D pattern to CSV for visualization. */
void write_pattern_csv(const std::string &path, double phi_deg,
                       const std::vector<FFPoint2D> &pat);

} // namespace vectorem

