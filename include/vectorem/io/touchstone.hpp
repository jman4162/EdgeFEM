#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include "vectorem/ports/port_eigensolve.hpp"

namespace vectorem {

/// Output format for Touchstone files.
enum class TouchstoneFormat {
  RI,  // Real/Imaginary
  MA,  // Magnitude/Angle (degrees)
  DB   // dB/Angle (degrees)
};

/// Options for Touchstone file export.
struct TouchstoneOptions {
  TouchstoneFormat format = TouchstoneFormat::RI;
  double z0 = 50.0;  // Reference impedance (ohms)
};

/// Write 2-port S-parameters to Touchstone file (legacy API).
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
void write_touchstone(const std::string &path, const std::vector<double> &freq,
                      const std::vector<SParams2> &data);
#pragma GCC diagnostic pop

/// Write N-port S-parameters to Touchstone file (.sNp format).
/// @param path Output file path (extension determined automatically)
/// @param freq Vector of frequencies in Hz
/// @param S_matrices Vector of NxN S-parameter matrices, one per frequency
/// @param opts Export options (format, reference impedance)
void write_touchstone_nport(const std::string &path,
                            const std::vector<double> &freq,
                            const std::vector<Eigen::MatrixXcd> &S_matrices,
                            const TouchstoneOptions &opts = TouchstoneOptions());

/// Get the appropriate Touchstone file extension for a given port count.
/// @param num_ports Number of ports (1-99)
/// @return Extension string (e.g., ".s2p", ".s3p", ".s4p")
std::string touchstone_extension(int num_ports);

} // namespace vectorem
