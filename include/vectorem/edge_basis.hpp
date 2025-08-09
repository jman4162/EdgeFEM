#pragma once

#include <array>
#include <Eigen/Core>

namespace vectorem {

Eigen::Matrix<double, 6, 3>
whitney_edge_curls(const std::array<Eigen::Vector3d, 4> &v);

Eigen::Matrix<double, 6, 6>
whitney_curl_curl_matrix(const std::array<Eigen::Vector3d, 4> &v);

Eigen::Matrix<double, 6, 6>
whitney_mass_matrix(const std::array<Eigen::Vector3d, 4> &v);

}
