#include "test_check.hpp"
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "edgefem/io/touchstone.hpp"

using namespace edgefem;

void test_touchstone_extension() {
  std::cout << "test_touchstone_extension..." << std::endl;

  CHECK(touchstone_extension(1) == ".s1p");
  CHECK(touchstone_extension(2) == ".s2p");
  CHECK(touchstone_extension(3) == ".s3p");
  CHECK(touchstone_extension(4) == ".s4p");
  CHECK(touchstone_extension(10) == ".s10p");

  std::cout << "  Extension for 4 ports: " << touchstone_extension(4)
            << std::endl;
}

void test_write_3port() {
  std::cout << "test_write_3port..." << std::endl;

  std::vector<double> freq = {1e9, 2e9, 3e9};
  std::vector<Eigen::MatrixXcd> S_matrices;

  for (size_t i = 0; i < freq.size(); ++i) {
    Eigen::MatrixXcd S(3, 3);
    // Create a simple passive network
    S << std::complex<double>(0.1, 0.0), std::complex<double>(-0.3, 0.1),
        std::complex<double>(-0.2, 0.0), std::complex<double>(-0.3, 0.1),
        std::complex<double>(0.1, 0.0), std::complex<double>(-0.3, 0.1),
        std::complex<double>(-0.2, 0.0), std::complex<double>(-0.3, 0.1),
        std::complex<double>(0.1, 0.0);
    S_matrices.push_back(S);
  }

  TouchstoneOptions opts;
  opts.format = TouchstoneFormat::RI;
  opts.z0 = 50.0;

  std::string path = "/tmp/test_3port.s3p";
  write_touchstone_nport(path, freq, S_matrices, opts);

  // Verify file was created
  std::ifstream ifs(path);
  CHECK(ifs.good());

  // Read first few lines
  std::string line;
  std::getline(ifs, line); // Comment
  std::getline(ifs, line); // Comment
  std::getline(ifs, line); // Option line

  CHECK(line.find("Hz S RI R 50") != std::string::npos);

  std::cout << "  Wrote 3-port file: " << path << std::endl;
}

void test_write_4port() {
  std::cout << "test_write_4port..." << std::endl;

  std::vector<double> freq = {1e9, 2e9};
  std::vector<Eigen::MatrixXcd> S_matrices;

  for (size_t i = 0; i < freq.size(); ++i) {
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(4, 4);
    // Diagonal reflection
    for (int j = 0; j < 4; ++j) {
      S(j, j) = std::complex<double>(0.1, 0.0);
    }
    // Transmission
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        if (j != k) {
          S(j, k) = std::complex<double>(-0.2, 0.05);
        }
      }
    }
    S_matrices.push_back(S);
  }

  TouchstoneOptions opts;
  opts.format = TouchstoneFormat::MA; // Magnitude/Angle format
  opts.z0 = 50.0;

  std::string path = "/tmp/test_4port.s4p";
  write_touchstone_nport(path, freq, S_matrices, opts);

  // Verify file was created
  std::ifstream ifs(path);
  CHECK(ifs.good());

  std::string line;
  std::getline(ifs, line);
  std::getline(ifs, line);
  std::getline(ifs, line); // Option line

  CHECK(line.find("Hz S MA R 50") != std::string::npos);

  std::cout << "  Wrote 4-port file with MA format: " << path << std::endl;
}

void test_db_format() {
  std::cout << "test_db_format..." << std::endl;

  std::vector<double> freq = {1e9};
  Eigen::MatrixXcd S(2, 2);
  S << std::complex<double>(0.1, 0.0), std::complex<double>(-0.9, 0.1),
      std::complex<double>(-0.9, 0.1), std::complex<double>(0.1, 0.0);

  std::vector<Eigen::MatrixXcd> S_matrices = {S};

  TouchstoneOptions opts;
  opts.format = TouchstoneFormat::DB;
  opts.z0 = 50.0;

  std::string path = "/tmp/test_2port_db.s2p";
  write_touchstone_nport(path, freq, S_matrices, opts);

  std::ifstream ifs(path);
  CHECK(ifs.good());

  std::string line;
  std::getline(ifs, line);
  std::getline(ifs, line);
  std::getline(ifs, line);

  CHECK(line.find("Hz S DB R 50") != std::string::npos);

  std::cout << "  Wrote 2-port file with dB format: " << path << std::endl;
}

void test_empty_data_throws() {
  std::cout << "test_empty_data_throws..." << std::endl;

  std::vector<double> freq;
  std::vector<Eigen::MatrixXcd> S_matrices;

  bool threw = false;
  try {
    write_touchstone_nport("/tmp/empty.s2p", freq, S_matrices);
  } catch (const std::runtime_error &) {
    threw = true;
  }

  CHECK(threw);
  std::cout << "  Empty data correctly threw exception" << std::endl;
}

void test_dimension_mismatch_throws() {
  std::cout << "test_dimension_mismatch_throws..." << std::endl;

  std::vector<double> freq = {1e9, 2e9, 3e9};
  Eigen::MatrixXcd S(2, 2);
  S << 0.1, -0.9, -0.9, 0.1;

  // Only one matrix but 3 frequencies
  std::vector<Eigen::MatrixXcd> S_matrices = {S};

  bool threw = false;
  try {
    write_touchstone_nport("/tmp/mismatch.s2p", freq, S_matrices);
  } catch (const std::runtime_error &) {
    threw = true;
  }

  CHECK(threw);
  std::cout << "  Dimension mismatch correctly threw exception" << std::endl;
}

// Data lines (option line and comments dropped), each split into numbers.
std::vector<std::vector<double>> data_lines(const std::string &path) {
  std::ifstream ifs(path);
  std::vector<std::vector<double>> out;
  std::string line;
  while (std::getline(ifs, line)) {
    const auto bang = line.find('!');
    if (bang != std::string::npos)
      line.resize(bang);
    if (line.find('#') != std::string::npos)
      continue;
    std::istringstream is(line);
    std::vector<double> nums;
    double v;
    while (is >> v)
      nums.push_back(v);
    if (!nums.empty())
      out.push_back(nums);
  }
  return out;
}

Eigen::MatrixXcd numbered(int n) {
  // S(i, j) = (10 (i + 1) + (j + 1)) + j 0.5: every entry distinct, so the
  // order on the line identifies it.
  Eigen::MatrixXcd S(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      S(i, j) = std::complex<double>(10.0 * (i + 1) + (j + 1), 0.5);
  return S;
}

// Touchstone 1.1: a 2-port line is S11 S21 S12 S22 (column order); for
// N >= 3 each matrix row starts on a new line, S(i,1) .. S(i,N), at most
// four pairs per line.
void test_data_layout() {
  std::cout << "Testing Touchstone data layout..." << std::endl;
  {
    const std::string path = "/tmp/edgefem_layout.s2p";
    write_touchstone_nport(path, {1e9}, {numbered(2)});
    auto lines = data_lines(path);
    CHECK(lines.size() == 1);
    if (lines.size() == 1) {
      CHECK(lines[0].size() == 9);
      const double want[4] = {11, 21, 12, 22};
      for (int k = 0; k < 4 && lines[0].size() == 9; ++k)
        CHECK(lines[0][1 + 2 * k] == want[k]);
    }
  }
  for (int n : {3, 5}) {
    const std::string path = "/tmp/edgefem_layout" + touchstone_extension(n);
    write_touchstone_nport(path, {1e9}, {numbered(n)});
    auto lines = data_lines(path);
    // Lines per row: ceil(n / 4); the first line also carries the frequency.
    const int per_row = (n + 3) / 4;
    CHECK(static_cast<int>(lines.size()) == n * per_row);
    size_t li = 0;
    for (int i = 0; i < n && li < lines.size(); ++i) {
      std::vector<double> row;
      for (int k = 0; k < per_row && li < lines.size(); ++k, ++li) {
        const auto &ln = lines[li];
        const size_t first = (i == 0 && k == 0) ? 1 : 0; // skip frequency
        for (size_t t = first; t < ln.size(); ++t)
          row.push_back(ln[t]);
      }
      CHECK(static_cast<int>(row.size()) == 2 * n);
      for (int j = 0; j < n && static_cast<int>(row.size()) == 2 * n; ++j)
        CHECK(row[2 * j] == 10.0 * (i + 1) + (j + 1));
    }
  }
}

int main() {
  test_touchstone_extension();
  test_write_3port();
  test_write_4port();
  test_db_format();
  test_empty_data_throws();
  test_dimension_mismatch_throws();
  test_data_layout();
  std::cout << "All Touchstone N-port tests passed!" << std::endl;
  return edgefem_test::exit_code();
}
