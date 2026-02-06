#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

#include "edgefem/io/touchstone.hpp"

using namespace edgefem;

void test_touchstone_extension() {
  std::cout << "test_touchstone_extension..." << std::endl;

  assert(touchstone_extension(1) == ".s1p");
  assert(touchstone_extension(2) == ".s2p");
  assert(touchstone_extension(3) == ".s3p");
  assert(touchstone_extension(4) == ".s4p");
  assert(touchstone_extension(10) == ".s10p");

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
  assert(ifs.good());

  // Read first few lines
  std::string line;
  std::getline(ifs, line); // Comment
  std::getline(ifs, line); // Comment
  std::getline(ifs, line); // Option line

  assert(line.find("Hz S RI R 50") != std::string::npos);

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
  assert(ifs.good());

  std::string line;
  std::getline(ifs, line);
  std::getline(ifs, line);
  std::getline(ifs, line); // Option line

  assert(line.find("Hz S MA R 50") != std::string::npos);

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
  assert(ifs.good());

  std::string line;
  std::getline(ifs, line);
  std::getline(ifs, line);
  std::getline(ifs, line);

  assert(line.find("Hz S DB R 50") != std::string::npos);

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

  assert(threw);
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

  assert(threw);
  std::cout << "  Dimension mismatch correctly threw exception" << std::endl;
}

int main() {
  test_touchstone_extension();
  test_write_3port();
  test_write_4port();
  test_db_format();
  test_empty_data_throws();
  test_dimension_mismatch_throws();
  std::cout << "All Touchstone N-port tests passed!" << std::endl;
  return 0;
}
