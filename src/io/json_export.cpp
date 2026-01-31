#include "vectorem/io/json_export.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace vectorem {

namespace {

std::string complex_to_json(std::complex<double> val) {
  std::ostringstream ss;
  ss << std::setprecision(12);
  ss << "[" << val.real() << ", " << val.imag() << "]";
  return ss.str();
}

std::string matrix_to_json(const Eigen::MatrixXcd &M, int indent = 4) {
  std::ostringstream ss;
  std::string pad(indent, ' ');
  std::string inner_pad(indent + 2, ' ');

  ss << "[\n";
  for (int i = 0; i < M.rows(); ++i) {
    ss << inner_pad << "[";
    for (int j = 0; j < M.cols(); ++j) {
      ss << complex_to_json(M(i, j));
      if (j < M.cols() - 1)
        ss << ", ";
    }
    ss << "]";
    if (i < M.rows() - 1)
      ss << ",";
    ss << "\n";
  }
  ss << pad << "]";
  return ss.str();
}

} // namespace

void write_coupling_json(const std::string &path,
                         const std::vector<CouplingMatrix> &data) {
  if (data.empty()) {
    throw std::runtime_error("Empty coupling data");
  }

  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  ofs << std::setprecision(12);
  ofs << "{\n";
  ofs << "  \"type\": \"coupling_matrix\",\n";
  ofs << "  \"version\": \"1.0\",\n";
  ofs << "  \"num_ports\": " << data[0].num_ports << ",\n";
  ofs << "  \"z0\": " << data[0].z0 << ",\n";
  ofs << "  \"num_frequencies\": " << data.size() << ",\n";
  ofs << "  \"data\": [\n";

  for (size_t i = 0; i < data.size(); ++i) {
    const auto &cm = data[i];
    ofs << "    {\n";
    ofs << "      \"frequency\": " << cm.frequency << ",\n";
    ofs << "      \"S\": " << matrix_to_json(cm.S, 6) << ",\n";
    ofs << "      \"Z\": " << matrix_to_json(cm.Z, 6) << ",\n";
    ofs << "      \"Y\": " << matrix_to_json(cm.Y, 6) << "\n";
    ofs << "    }";
    if (i < data.size() - 1)
      ofs << ",";
    ofs << "\n";
  }

  ofs << "  ]\n";
  ofs << "}\n";
}

void write_sparams_json(const std::string &path,
                        const std::vector<double> &freq,
                        const std::vector<Eigen::MatrixXcd> &S) {
  if (freq.empty() || S.empty()) {
    throw std::runtime_error("Empty frequency or S-parameter data");
  }

  if (freq.size() != S.size()) {
    throw std::runtime_error("Frequency and S-matrix count mismatch");
  }

  std::ofstream ofs(path);
  if (!ofs) {
    throw std::runtime_error("Cannot open file for writing: " + path);
  }

  int num_ports = static_cast<int>(S[0].rows());

  ofs << std::setprecision(12);
  ofs << "{\n";
  ofs << "  \"type\": \"s_parameters\",\n";
  ofs << "  \"version\": \"1.0\",\n";
  ofs << "  \"num_ports\": " << num_ports << ",\n";
  ofs << "  \"num_frequencies\": " << freq.size() << ",\n";
  ofs << "  \"frequencies\": [";
  for (size_t i = 0; i < freq.size(); ++i) {
    ofs << freq[i];
    if (i < freq.size() - 1)
      ofs << ", ";
  }
  ofs << "],\n";
  ofs << "  \"S_matrices\": [\n";

  for (size_t i = 0; i < S.size(); ++i) {
    ofs << "    " << matrix_to_json(S[i], 4);
    if (i < S.size() - 1)
      ofs << ",";
    ofs << "\n";
  }

  ofs << "  ]\n";
  ofs << "}\n";
}

} // namespace vectorem
