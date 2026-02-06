#include "edgefem/array/coupling_extract.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace edgefem {

namespace {

double to_db(double linear) {
  return 20.0 * std::log10(std::max(linear, 1e-20));
}

double to_degrees(double radians) { return radians * 180.0 / M_PI; }

} // namespace

std::vector<CouplingCoefficient>
extract_all_couplings(const Eigen::MatrixXcd &S,
                      const std::vector<Eigen::Vector3d> &positions) {

  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }

  int n = static_cast<int>(S.rows());
  bool have_positions = (static_cast<int>(positions.size()) == n);

  std::vector<CouplingCoefficient> result;
  result.reserve(n * (n - 1));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;

      CouplingCoefficient coef;
      coef.port_i = i;
      coef.port_j = j;
      coef.S_ij = S(i, j);
      coef.magnitude_db = to_db(std::abs(S(i, j)));
      coef.phase_deg = to_degrees(std::arg(S(i, j)));

      if (have_positions) {
        coef.distance = (positions[i] - positions[j]).norm();
      } else {
        coef.distance = 0.0;
      }

      result.push_back(coef);
    }
  }

  return result;
}

std::map<int, double>
coupling_vs_separation(const Eigen::MatrixXcd &S,
                       const std::vector<Eigen::Vector3d> &positions) {

  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }

  int n = static_cast<int>(S.rows());
  if (static_cast<int>(positions.size()) != n) {
    throw std::runtime_error("Position count must match S-matrix dimension");
  }

  // Find min/max distances to set bin size
  double min_dist = std::numeric_limits<double>::max();
  double max_dist = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double d = (positions[i] - positions[j]).norm();
      min_dist = std::min(min_dist, d);
      max_dist = std::max(max_dist, d);
    }
  }

  if (max_dist < 1e-12) {
    return {};
  }

  // Bin distances (10 bins)
  double bin_size = (max_dist - min_dist) / 10.0;
  if (bin_size < 1e-12)
    bin_size = max_dist / 10.0;

  std::map<int, std::vector<double>> bins;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double d = (positions[i] - positions[j]).norm();
      int bin = static_cast<int>((d - min_dist) / bin_size);

      // Average of both directions (should be reciprocal for passive network)
      double coupling_db =
          0.5 * (to_db(std::abs(S(i, j))) + to_db(std::abs(S(j, i))));
      bins[bin].push_back(coupling_db);
    }
  }

  // Compute average per bin
  std::map<int, double> result;
  for (const auto &kv : bins) {
    double sum = 0.0;
    for (double v : kv.second) {
      sum += v;
    }
    result[kv.first] = sum / kv.second.size();
  }

  return result;
}

std::vector<CouplingCoefficient>
extract_nearest_neighbor_coupling(const Eigen::MatrixXcd &S,
                                  const std::vector<Eigen::Vector3d> &positions,
                                  int max_neighbors) {

  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }

  int n = static_cast<int>(S.rows());
  if (static_cast<int>(positions.size()) != n) {
    throw std::runtime_error("Position count must match S-matrix dimension");
  }

  std::vector<CouplingCoefficient> result;

  for (int i = 0; i < n; ++i) {
    // Find distances to all other elements
    std::vector<std::pair<double, int>> distances;
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;
      double d = (positions[i] - positions[j]).norm();
      distances.push_back({d, j});
    }

    // Sort by distance
    std::sort(distances.begin(), distances.end());

    // Take up to max_neighbors nearest
    int count = std::min(max_neighbors, static_cast<int>(distances.size()));
    for (int k = 0; k < count; ++k) {
      int j = distances[k].second;

      CouplingCoefficient coef;
      coef.port_i = i;
      coef.port_j = j;
      coef.S_ij = S(i, j);
      coef.magnitude_db = to_db(std::abs(S(i, j)));
      coef.phase_deg = to_degrees(std::arg(S(i, j)));
      coef.distance = distances[k].first;

      result.push_back(coef);
    }
  }

  return result;
}

CouplingStats compute_coupling_stats(const Eigen::MatrixXcd &S) {
  if (S.rows() != S.cols()) {
    throw std::runtime_error("S-matrix must be square");
  }

  int n = static_cast<int>(S.rows());
  CouplingStats stats;

  stats.max_coupling_db = -std::numeric_limits<double>::infinity();
  stats.min_coupling_db = std::numeric_limits<double>::infinity();
  stats.max_coupling_i = 0;
  stats.max_coupling_j = 0;

  double sum = 0.0;
  double sum_sq = 0.0;
  int count = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;

      double coupling_db = to_db(std::abs(S(i, j)));

      if (coupling_db > stats.max_coupling_db) {
        stats.max_coupling_db = coupling_db;
        stats.max_coupling_i = i;
        stats.max_coupling_j = j;
      }
      if (coupling_db < stats.min_coupling_db) {
        stats.min_coupling_db = coupling_db;
      }

      sum += coupling_db;
      sum_sq += coupling_db * coupling_db;
      count++;
    }
  }

  if (count > 0) {
    stats.mean_coupling_db = sum / count;
    double variance =
        (sum_sq / count) - (stats.mean_coupling_db * stats.mean_coupling_db);
    stats.std_coupling_db = std::sqrt(std::max(variance, 0.0));
  } else {
    stats.mean_coupling_db = 0.0;
    stats.std_coupling_db = 0.0;
  }

  return stats;
}

} // namespace edgefem
