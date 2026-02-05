#include "vectorem/post/ntf.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

namespace vectorem {

static constexpr double Z0 = 376.730313668; // free-space impedance [ohm]

// -----------------------------------------------------------------------------
// FFPattern3D member functions
// -----------------------------------------------------------------------------

Eigen::MatrixXd FFPattern3D::total_magnitude() const {
  Eigen::MatrixXd mag(E_theta.rows(), E_theta.cols());
  for (int i = 0; i < E_theta.rows(); ++i) {
    for (int j = 0; j < E_theta.cols(); ++j) {
      double eth2 = std::norm(E_theta(i, j));
      double eph2 = std::norm(E_phi(i, j));
      mag(i, j) = std::sqrt(eth2 + eph2);
    }
  }
  return mag;
}

Eigen::MatrixXd FFPattern3D::power_pattern() const {
  Eigen::MatrixXd pwr(E_theta.rows(), E_theta.cols());
  for (int i = 0; i < E_theta.rows(); ++i) {
    for (int j = 0; j < E_theta.cols(); ++j) {
      pwr(i, j) = std::norm(E_theta(i, j)) + std::norm(E_phi(i, j));
    }
  }
  return pwr;
}

Eigen::MatrixXd FFPattern3D::pattern_dB() const {
  Eigen::MatrixXd pwr = power_pattern();
  double max_pwr = pwr.maxCoeff();
  if (max_pwr < std::numeric_limits<double>::epsilon()) {
    return Eigen::MatrixXd::Constant(pwr.rows(), pwr.cols(), -200.0);
  }
  Eigen::MatrixXd dB(pwr.rows(), pwr.cols());
  for (int i = 0; i < pwr.rows(); ++i) {
    for (int j = 0; j < pwr.cols(); ++j) {
      double ratio = pwr(i, j) / max_pwr;
      dB(i, j) = (ratio > 1e-20) ? 10.0 * std::log10(ratio) : -200.0;
    }
  }
  return dB;
}

// -----------------------------------------------------------------------------
// 2D Stratton-Chu (existing)
// -----------------------------------------------------------------------------

std::vector<NTFPoint2D> stratton_chu_2d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    double phi_rad, double k0) {
  std::vector<NTFPoint2D> out;
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

    NTFPoint2D p;
    p.theta_deg = th * 180.0 / M_PI;
    p.e_theta = Efar.dot(th_hat);
    p.e_phi = Efar.dot(ph_hat);
    out.push_back(p);
  }
  return out;
}

// -----------------------------------------------------------------------------
// 3D Stratton-Chu
// -----------------------------------------------------------------------------

FFPattern3D stratton_chu_3d(
    const std::vector<Eigen::Vector3d> &r,
    const std::vector<Eigen::Vector3d> &n,
    const std::vector<Eigen::Vector3cd> &E,
    const std::vector<Eigen::Vector3cd> &H,
    const std::vector<double> &area,
    const std::vector<double> &theta_rad,
    const std::vector<double> &phi_rad,
    double k0) {
  const int Nth = static_cast<int>(theta_rad.size());
  const int Nph = static_cast<int>(phi_rad.size());

  FFPattern3D pattern;
  pattern.theta_grid.resize(Nth, Nph);
  pattern.phi_grid.resize(Nth, Nph);
  pattern.E_theta.resize(Nth, Nph);
  pattern.E_phi.resize(Nth, Nph);

  // Build grids
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      pattern.theta_grid(i, j) = theta_rad[i];
      pattern.phi_grid(i, j) = phi_rad[j];
    }
  }

  const std::complex<double> j_imag(0.0, 1.0);

  // Compute far-field at each (theta, phi) sample point
  for (int ti = 0; ti < Nth; ++ti) {
    double th = theta_rad[ti];
    double sin_th = std::sin(th);
    double cos_th = std::cos(th);

    for (int pi = 0; pi < Nph; ++pi) {
      double ph = phi_rad[pi];
      double sin_ph = std::sin(ph);
      double cos_ph = std::cos(ph);

      // Spherical unit vectors
      Eigen::Vector3d rhat(sin_th * cos_ph, sin_th * sin_ph, cos_th);
      Eigen::Vector3d th_hat(cos_th * cos_ph, cos_th * sin_ph, -sin_th);
      Eigen::Vector3d ph_hat(-sin_ph, cos_ph, 0.0);

      // Integrate over Huygens surface
      Eigen::Vector3cd Efar = Eigen::Vector3cd::Zero();
      for (size_t s = 0; s < r.size(); ++s) {
        // Equivalent surface currents
        Eigen::Vector3cd J = n[s].cross(H[s]);        // Electric current
        Eigen::Vector3cd M = -n[s].cross(E[s]);       // Magnetic current

        // Phase factor for this surface point
        std::complex<double> phase = std::exp(-j_imag * k0 * rhat.dot(r[s]));

        // Far-field contribution from this patch
        // E_far = jk0 * (rhat x M) x rhat - Z0 * rhat x J
        Eigen::Vector3cd term = j_imag * k0 * (rhat.cross(M)).cross(rhat) -
                                Z0 * rhat.cross(J);
        Efar += term * phase * area[s];
      }

      // Project onto spherical basis vectors
      pattern.E_theta(ti, pi) = Efar.dot(th_hat);
      pattern.E_phi(ti, pi) = Efar.dot(ph_hat);
    }
  }

  return pattern;
}

// -----------------------------------------------------------------------------
// Directivity and gain computations
// -----------------------------------------------------------------------------

double compute_directivity(const FFPattern3D &pattern) {
  // Directivity D = 4*pi * U_max / P_rad
  // U = r^2 * S = |E|^2 / (2*Z0) is radiation intensity
  // P_rad = integral of U over sphere

  const Eigen::MatrixXd pwr = pattern.power_pattern();  // |E_th|^2 + |E_ph|^2
  const double U_max = pwr.maxCoeff() / (2.0 * Z0);

  if (U_max < std::numeric_limits<double>::epsilon()) {
    return 0.0;
  }

  // Numerical integration over sphere using trapezoidal rule
  // P_rad = integral(U * sin(theta) dtheta dphi)
  const int Nth = static_cast<int>(pattern.theta_grid.rows());
  const int Nph = static_cast<int>(pattern.theta_grid.cols());

  if (Nth < 2 || Nph < 2) {
    return 1.0;  // Not enough points for integration
  }

  // Compute dtheta and dphi (assuming uniform spacing)
  double dtheta = (Nth > 1) ?
      (pattern.theta_grid(Nth - 1, 0) - pattern.theta_grid(0, 0)) / (Nth - 1) : M_PI;
  double dphi = (Nph > 1) ?
      (pattern.phi_grid(0, Nph - 1) - pattern.phi_grid(0, 0)) / (Nph - 1) : 2.0 * M_PI;

  double P_rad = 0.0;
  for (int i = 0; i < Nth; ++i) {
    double theta = pattern.theta_grid(i, 0);
    double sin_th = std::sin(theta);
    double wt_th = (i == 0 || i == Nth - 1) ? 0.5 : 1.0;  // Trapezoidal weight

    for (int j = 0; j < Nph; ++j) {
      double wt_ph = (j == 0 || j == Nph - 1) ? 0.5 : 1.0;
      double U = pwr(i, j) / (2.0 * Z0);
      P_rad += U * sin_th * wt_th * wt_ph;
    }
  }
  P_rad *= dtheta * dphi;

  if (P_rad < std::numeric_limits<double>::epsilon()) {
    return 0.0;
  }

  return 4.0 * M_PI * U_max / P_rad;
}

double compute_max_gain(const FFPattern3D &pattern, double efficiency) {
  return efficiency * compute_directivity(pattern);
}

std::pair<double, double> compute_hpbw(const FFPattern3D &pattern) {
  const Eigen::MatrixXd pwr = pattern.power_pattern();
  double max_pwr = pwr.maxCoeff();

  if (max_pwr < std::numeric_limits<double>::epsilon()) {
    return {0.0, 0.0};
  }

  double half_power = max_pwr / 2.0;
  const int Nth = static_cast<int>(pattern.theta_grid.rows());
  const int Nph = static_cast<int>(pattern.theta_grid.cols());

  // Find maximum location
  int max_ti = 0, max_pi = 0;
  pwr.maxCoeff(&max_ti, &max_pi);

  // E-plane: scan theta at phi = phi_max
  double e_plane_hpbw = 0.0;
  {
    int lower = max_ti, upper = max_ti;
    while (lower > 0 && pwr(lower, max_pi) > half_power) --lower;
    while (upper < Nth - 1 && pwr(upper, max_pi) > half_power) ++upper;

    double theta_low = pattern.theta_grid(lower, max_pi);
    double theta_high = pattern.theta_grid(upper, max_pi);
    e_plane_hpbw = (theta_high - theta_low) * 180.0 / M_PI;
  }

  // H-plane: scan phi at theta = theta_max
  double h_plane_hpbw = 0.0;
  {
    int lower = max_pi, upper = max_pi;
    while (lower > 0 && pwr(max_ti, lower) > half_power) --lower;
    while (upper < Nph - 1 && pwr(max_ti, upper) > half_power) ++upper;

    double phi_low = pattern.phi_grid(max_ti, lower);
    double phi_high = pattern.phi_grid(max_ti, upper);
    // Account for sin(theta) scaling in H-plane
    double sin_th = std::sin(pattern.theta_grid(max_ti, max_pi));
    h_plane_hpbw = (phi_high - phi_low) * sin_th * 180.0 / M_PI;
  }

  return {e_plane_hpbw, h_plane_hpbw};
}

// -----------------------------------------------------------------------------
// CSV/VTK export
// -----------------------------------------------------------------------------

void write_pattern_csv(const std::string &path, double phi_deg,
                       const std::vector<NTFPoint2D> &pat) {
  std::ofstream f(path);
  f << "theta_deg,phi_deg,eth_mag,eph_mag\n";
  for (const auto &p : pat) {
    double eth = std::abs(p.e_theta);
    double eph = std::abs(p.e_phi);
    f << p.theta_deg << ',' << phi_deg << ',' << eth << ',' << eph << "\n";
  }
}

void write_pattern_3d_csv(const std::string &path, const FFPattern3D &pattern) {
  std::ofstream f(path);
  f << "theta_deg,phi_deg,E_theta_mag,E_theta_phase_deg,E_phi_mag,E_phi_phase_deg,total_dB\n";

  const Eigen::MatrixXd dB = pattern.pattern_dB();

  for (int i = 0; i < pattern.theta_grid.rows(); ++i) {
    for (int j = 0; j < pattern.phi_grid.cols(); ++j) {
      double th_deg = pattern.theta_grid(i, j) * 180.0 / M_PI;
      double ph_deg = pattern.phi_grid(i, j) * 180.0 / M_PI;

      std::complex<double> eth = pattern.E_theta(i, j);
      std::complex<double> eph = pattern.E_phi(i, j);

      f << th_deg << ','
        << ph_deg << ','
        << std::abs(eth) << ','
        << std::arg(eth) * 180.0 / M_PI << ','
        << std::abs(eph) << ','
        << std::arg(eph) * 180.0 / M_PI << ','
        << dB(i, j) << '\n';
    }
  }
}

void write_pattern_3d_vtk(const std::string &path, const FFPattern3D &pattern,
                          double radius) {
  std::ofstream f(path);

  const int Nth = static_cast<int>(pattern.theta_grid.rows());
  const int Nph = static_cast<int>(pattern.phi_grid.cols());
  const int n_points = Nth * Nph;
  const int n_cells = (Nth - 1) * (Nph - 1);

  // VTK XML Unstructured Grid header
  f << "<?xml version=\"1.0\"?>\n";
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  f << "  <UnstructuredGrid>\n";
  f << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";

  // Points (spherical surface)
  f << "      <Points>\n";
  f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      double th = pattern.theta_grid(i, j);
      double ph = pattern.phi_grid(i, j);
      double x = radius * std::sin(th) * std::cos(ph);
      double y = radius * std::sin(th) * std::sin(ph);
      double z = radius * std::cos(th);
      f << "          " << x << " " << y << " " << z << "\n";
    }
  }
  f << "        </DataArray>\n";
  f << "      </Points>\n";

  // Cells (quads)
  f << "      <Cells>\n";
  f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int i = 0; i < Nth - 1; ++i) {
    for (int j = 0; j < Nph - 1; ++j) {
      int p0 = i * Nph + j;
      int p1 = i * Nph + (j + 1);
      int p2 = (i + 1) * Nph + (j + 1);
      int p3 = (i + 1) * Nph + j;
      f << "          " << p0 << " " << p1 << " " << p2 << " " << p3 << "\n";
    }
  }
  f << "        </DataArray>\n";
  f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int c = 1; c <= n_cells; ++c) {
    f << "          " << (c * 4) << "\n";
  }
  f << "        </DataArray>\n";
  f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int c = 0; c < n_cells; ++c) {
    f << "          9\n";  // VTK_QUAD = 9
  }
  f << "        </DataArray>\n";
  f << "      </Cells>\n";

  // Point data
  const Eigen::MatrixXd dB = pattern.pattern_dB();
  const Eigen::MatrixXd mag = pattern.total_magnitude();

  f << "      <PointData Scalars=\"Pattern_dB\">\n";

  // Pattern in dB
  f << "        <DataArray type=\"Float64\" Name=\"Pattern_dB\" format=\"ascii\">\n";
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      f << "          " << dB(i, j) << "\n";
    }
  }
  f << "        </DataArray>\n";

  // Total magnitude
  f << "        <DataArray type=\"Float64\" Name=\"E_magnitude\" format=\"ascii\">\n";
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      f << "          " << mag(i, j) << "\n";
    }
  }
  f << "        </DataArray>\n";

  // E_theta magnitude
  f << "        <DataArray type=\"Float64\" Name=\"E_theta_mag\" format=\"ascii\">\n";
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      f << "          " << std::abs(pattern.E_theta(i, j)) << "\n";
    }
  }
  f << "        </DataArray>\n";

  // E_phi magnitude
  f << "        <DataArray type=\"Float64\" Name=\"E_phi_mag\" format=\"ascii\">\n";
  for (int i = 0; i < Nth; ++i) {
    for (int j = 0; j < Nph; ++j) {
      f << "          " << std::abs(pattern.E_phi(i, j)) << "\n";
    }
  }
  f << "        </DataArray>\n";

  f << "      </PointData>\n";
  f << "    </Piece>\n";
  f << "  </UnstructuredGrid>\n";
  f << "</VTKFile>\n";
}

} // namespace vectorem

