#include "vectorem/io/array_export.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "vectorem/array/coupling_extract.hpp"

namespace vectorem {

namespace {

std::string to_json_array(const Eigen::Vector3d &v) {
  std::ostringstream ss;
  ss << std::setprecision(12);
  ss << "[" << v.x() << ", " << v.y() << ", " << v.z() << "]";
  return ss.str();
}

std::string complex_to_json(std::complex<double> val) {
  std::ostringstream ss;
  ss << std::setprecision(12);
  ss << "[" << val.real() << ", " << val.imag() << "]";
  return ss.str();
}

double to_degrees(double radians) {
  return radians * 180.0 / M_PI;
}

} // namespace

void export_array_package_json(const std::string &path, const ArrayPackage &pkg) {
  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  ofs << std::setprecision(12);
  ofs << "{\n";
  ofs << "  \"format\": \"VectorEM_ArrayPackage\",\n";
  ofs << "  \"version\": \"1.0\",\n";
  ofs << "  \"model_name\": \"" << pkg.model_name << "\",\n";
  ofs << "  \"reference_impedance\": " << pkg.reference_impedance << ",\n";

  // Element positions
  ofs << "  \"num_elements\": " << pkg.element_positions.size() << ",\n";
  ofs << "  \"element_positions\": [\n";
  for (size_t i = 0; i < pkg.element_positions.size(); ++i) {
    ofs << "    " << to_json_array(pkg.element_positions[i]);
    if (i < pkg.element_positions.size() - 1) ofs << ",";
    ofs << "\n";
  }
  ofs << "  ],\n";

  // Frequencies
  ofs << "  \"num_frequencies\": " << pkg.frequencies.size() << ",\n";
  ofs << "  \"frequencies\": [";
  for (size_t i = 0; i < pkg.frequencies.size(); ++i) {
    ofs << pkg.frequencies[i];
    if (i < pkg.frequencies.size() - 1) ofs << ", ";
  }
  ofs << "],\n";

  // S-matrices
  ofs << "  \"S_matrices\": [\n";
  for (size_t fi = 0; fi < pkg.S_matrices.size(); ++fi) {
    const auto &S = pkg.S_matrices[fi];
    ofs << "    {\n";
    ofs << "      \"frequency\": " << pkg.frequencies[fi] << ",\n";
    ofs << "      \"data\": [\n";
    for (int i = 0; i < S.rows(); ++i) {
      ofs << "        [";
      for (int j = 0; j < S.cols(); ++j) {
        ofs << complex_to_json(S(i, j));
        if (j < S.cols() - 1) ofs << ", ";
      }
      ofs << "]";
      if (i < S.rows() - 1) ofs << ",";
      ofs << "\n";
    }
    ofs << "      ]\n";
    ofs << "    }";
    if (fi < pkg.S_matrices.size() - 1) ofs << ",";
    ofs << "\n";
  }
  ofs << "  ],\n";

  // Embedded patterns (simplified - just summary)
  ofs << "  \"has_embedded_patterns\": " << (pkg.embedded_patterns.empty() ? "false" : "true") << ",\n";
  ofs << "  \"num_pattern_frequencies\": " << pkg.embedded_patterns.size() << ",\n";

  // Scan data
  ofs << "  \"has_scan_data\": " << (pkg.scan_data.empty() ? "false" : "true") << ",\n";
  ofs << "  \"num_scan_points\": " << pkg.scan_data.size();

  if (!pkg.scan_data.empty()) {
    ofs << ",\n";
    ofs << "  \"scan_angles\": [\n";
    for (size_t i = 0; i < pkg.scan_data.size(); ++i) {
      ofs << "    {\"theta\": " << to_degrees(pkg.scan_data[i].theta)
          << ", \"phi\": " << to_degrees(pkg.scan_data[i].phi)
          << ", \"scan_loss_db\": " << pkg.scan_data[i].scan_loss_db << "}";
      if (i < pkg.scan_data.size() - 1) ofs << ",";
      ofs << "\n";
    }
    ofs << "  ]\n";
  } else {
    ofs << "\n";
  }

  ofs << "}\n";
}

void export_patterns_csv(const std::string &path,
                         const std::vector<EmbeddedPattern> &patterns) {
  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  ofs << std::setprecision(12);
  ofs << "port_index,theta_deg,phi0_mag,phi0_phase_deg,phi90_mag,phi90_phase_deg\n";

  for (const auto &pattern : patterns) {
    size_t n = std::min(pattern.pattern_phi0.size(), pattern.pattern_phi90.size());
    for (size_t i = 0; i < n; ++i) {
      const auto &p0 = pattern.pattern_phi0[i];
      const auto &p90 = pattern.pattern_phi90[i];
      ofs << pattern.port_index << ","
          << to_degrees(p0.theta) << ","
          << p0.magnitude << ","
          << to_degrees(p0.phase) << ","
          << p90.magnitude << ","
          << to_degrees(p90.phase) << "\n";
    }
  }
}

void export_coupling_csv(const std::string &path,
                         const Eigen::MatrixXcd &S,
                         const std::vector<Eigen::Vector3d> &positions) {
  auto couplings = extract_all_couplings(S, positions);

  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  ofs << std::setprecision(12);
  ofs << "port_i,port_j,S_real,S_imag,magnitude_db,phase_deg,distance_m\n";

  for (const auto &c : couplings) {
    ofs << c.port_i << ","
        << c.port_j << ","
        << c.S_ij.real() << ","
        << c.S_ij.imag() << ","
        << c.magnitude_db << ","
        << c.phase_deg << ","
        << c.distance << "\n";
  }
}

void export_scan_csv(const std::string &path,
                     const std::vector<ActiveImpedanceResult> &scan_results) {
  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  ofs << std::setprecision(12);
  ofs << "theta_deg,phi_deg,element_idx,Gamma_real,Gamma_imag,Z_real,Z_imag,VSWR\n";

  for (const auto &result : scan_results) {
    int n = static_cast<int>(result.Gamma_active.size());
    for (int i = 0; i < n; ++i) {
      double vswr = active_vswr(result.Gamma_active(i));
      ofs << to_degrees(result.theta) << ","
          << to_degrees(result.phi) << ","
          << i << ","
          << result.Gamma_active(i).real() << ","
          << result.Gamma_active(i).imag() << ","
          << result.Z_active(i).real() << ","
          << result.Z_active(i).imag() << ","
          << vswr << "\n";
    }
  }
}

} // namespace vectorem
