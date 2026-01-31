#include "vectorem/post/ntf.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <fstream>

namespace vectorem {

static constexpr double Z0 = 376.730313668; // free-space impedance [ohm]

std::vector<FFPoint2D> stratton_chu_2d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    double phi_rad, double k0) {
  std::vector<FFPoint2D> out;
  out.reserve(theta_rad.size());

  for (double th : theta_rad) {
    Eigen::Vector3d rhat(std::sin(th) * std::cos(phi_rad),
                         std::sin(th) * std::sin(phi_rad), std::cos(th));
    Eigen::Vector3d th_hat(std::cos(th) * std::cos(phi_rad),
                           std::cos(th) * std::sin(phi_rad), -std::sin(th));
    Eigen::Vector3d ph_hat(-std::sin(phi_rad), std::cos(phi_rad), 0.0);

    Eigen::Vector3cd Efar = Eigen::Vector3cd::Zero();
    const std::complex<double> j(0.0, 1.0);

    for (size_t i = 0; i < r.size(); ++i) {
      Eigen::Vector3cd J = n[i].cross(H[i]);
      Eigen::Vector3cd M = -n[i].cross(E[i]);
      std::complex<double> phase = std::exp(-j * k0 * rhat.dot(r[i]));
      Eigen::Vector3cd term = j * k0 * (rhat.cross(M)).cross(rhat) -
                              Z0 * rhat.cross(J);
      Efar += term * phase * area[i];
    }

    FFPoint2D p;
    p.theta_deg = th * 180.0 / M_PI;
    p.e_theta = Efar.dot(th_hat);
    p.e_phi = Efar.dot(ph_hat);
    out.push_back(p);
  }
  return out;
}

void write_pattern_csv(const std::string &path, double phi_deg,
                       const std::vector<FFPoint2D> &pat) {
  std::ofstream f(path);
  f << "theta_deg,phi_deg,eth_mag,eph_mag\n";
  for (const auto &p : pat) {
    double eth = std::abs(p.e_theta);
    double eph = std::abs(p.e_phi);
    f << p.theta_deg << ',' << phi_deg << ',' << eth << ',' << eph << "\n";
  }
}

} // namespace vectorem

