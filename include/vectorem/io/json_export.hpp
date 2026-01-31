#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include "vectorem/coupling.hpp"

namespace vectorem {

/// Write coupling matrix data to JSON file.
/// @param path Output file path
/// @param data Vector of CouplingMatrix objects (one per frequency)
void write_coupling_json(const std::string &path,
                         const std::vector<CouplingMatrix> &data);

/// Write S-parameter data to JSON file.
/// @param path Output file path
/// @param freq Vector of frequencies in Hz
/// @param S Vector of S-parameter matrices (one per frequency)
void write_sparams_json(const std::string &path,
                        const std::vector<double> &freq,
                        const std::vector<Eigen::MatrixXcd> &S);

} // namespace vectorem
