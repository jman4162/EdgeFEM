#include <cassert>
#include <vector>

#include "vectorem/post/ntf.hpp"

using namespace vectorem;

int main() {
  std::vector<Eigen::Vector3d> r = {Eigen::Vector3d::Zero()};
  std::vector<Eigen::Vector3d> n = {Eigen::Vector3d::UnitZ()};
  std::vector<Eigen::Vector3cd> E = {Eigen::Vector3cd::Zero()};
  std::vector<Eigen::Vector3cd> H = {Eigen::Vector3cd::Zero()};
  std::vector<double> theta = {0.0, M_PI / 2};
  auto pat = stratton_chu_2d(r, n, E, H, theta, 0.0, 2 * M_PI);
  assert(pat.size() == 2);
  assert(std::abs(pat[0].e_theta) < 1e-12);
  return 0;
}
