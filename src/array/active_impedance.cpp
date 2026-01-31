#include "vectorem/array/active_impedance.hpp"

#include <cmath>
#include <stdexcept>

namespace vectorem {

Eigen::VectorXcd uniform_scan_excitation(int num_elements,
                                         const std::vector<Eigen::Vector3d> &positions,
                                         double theta, double phi, double k0) {
  if (static_cast<int>(positions.size()) != num_elements) {
    throw std::runtime_error("Position vector size must match num_elements");
  }

  Eigen::VectorXcd excitation(num_elements);

  // Scan direction unit vector
  Eigen::Vector3d k_hat;
  k_hat.x() = std::sin(theta) * std::cos(phi);
  k_hat.y() = std::sin(theta) * std::sin(phi);
  k_hat.z() = std::cos(theta);

  // Phase for each element: exp(-j * k0 * k_hat · position)
  for (int i = 0; i < num_elements; ++i) {
    double phase = -k0 * k_hat.dot(positions[i]);
    excitation(i) = std::complex<double>(std::cos(phase), std::sin(phase));
  }

  return excitation;
}

ActiveImpedanceResult compute_active_impedance(const CouplingMatrix &coupling,
                                                const Eigen::VectorXcd &excitation) {
  if (coupling.S.rows() != excitation.size()) {
    throw std::runtime_error("Excitation size must match S-matrix dimension");
  }

  ActiveImpedanceResult result;
  result.theta = 0.0;  // Caller should set these
  result.phi = 0.0;

  int n = coupling.num_ports;

  // Compute active reflection coefficients
  result.Gamma_active = active_reflection(coupling.S, excitation);

  // Compute active impedances: Z_active = Z0 * (1 + Gamma) / (1 - Gamma)
  result.Z_active.resize(n);
  for (int i = 0; i < n; ++i) {
    std::complex<double> gamma = result.Gamma_active(i);
    std::complex<double> denom = 1.0 - gamma;

    if (std::abs(denom) < 1e-12) {
      // Near total reflection
      result.Z_active(i) = std::complex<double>(1e12, 0);
    } else {
      result.Z_active(i) = coupling.z0 * (1.0 + gamma) / denom;
    }
  }

  // Compute scan loss (average power loss relative to matched condition)
  double total_reflected_power = 0.0;
  double total_incident_power = 0.0;
  for (int i = 0; i < n; ++i) {
    double incident = std::norm(excitation(i));
    double reflected = std::norm(result.Gamma_active(i)) * incident;
    total_incident_power += incident;
    total_reflected_power += reflected;
  }

  if (total_incident_power > 1e-20) {
    double mismatch_efficiency = 1.0 - (total_reflected_power / total_incident_power);
    result.scan_loss_db = -10.0 * std::log10(std::max(mismatch_efficiency, 1e-20));
  } else {
    result.scan_loss_db = 0.0;
  }

  return result;
}

std::vector<ActiveImpedanceResult> active_impedance_scan(
    const CouplingMatrix &coupling,
    const std::vector<Eigen::Vector3d> &positions,
    const std::vector<double> &theta_rad,
    const std::vector<double> &phi_rad,
    double k0) {

  int num_elements = coupling.num_ports;

  if (static_cast<int>(positions.size()) != num_elements) {
    throw std::runtime_error("Position count must match number of ports");
  }

  std::vector<ActiveImpedanceResult> results;
  results.reserve(theta_rad.size() * phi_rad.size());

  for (double theta : theta_rad) {
    for (double phi : phi_rad) {
      Eigen::VectorXcd excitation = uniform_scan_excitation(
          num_elements, positions, theta, phi, k0);

      ActiveImpedanceResult result = compute_active_impedance(coupling, excitation);
      result.theta = theta;
      result.phi = phi;

      results.push_back(result);
    }
  }

  return results;
}

double active_vswr(std::complex<double> Gamma_active) {
  double rho = std::abs(Gamma_active);
  if (rho >= 1.0 - 1e-12) {
    return 1e12;  // Effectively infinite VSWR
  }
  return (1.0 + rho) / (1.0 - rho);
}

std::map<int, std::vector<double>> find_scan_blindness(
    const std::vector<ActiveImpedanceResult> &scan_results,
    double vswr_threshold) {

  std::map<int, std::vector<double>> blindness;

  for (const auto &result : scan_results) {
    int n = static_cast<int>(result.Gamma_active.size());
    for (int i = 0; i < n; ++i) {
      double vswr = active_vswr(result.Gamma_active(i));
      if (vswr > vswr_threshold) {
        blindness[i].push_back(result.theta);
      }
    }
  }

  return blindness;
}

} // namespace vectorem
