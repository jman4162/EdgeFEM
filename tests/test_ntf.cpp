#include <cassert>
#include <vector>
#include <cmath>

#include "vectorem/post/ntf.hpp"

using namespace vectorem;

int main() {
  // Zero-current sanity check
  {
    std::vector<Eigen::Vector3d> r = {Eigen::Vector3d::Zero()};
    std::vector<Eigen::Vector3d> n = {Eigen::Vector3d::UnitZ()};
    std::vector<Eigen::Vector3cd> E = {Eigen::Vector3cd::Zero()};
    std::vector<Eigen::Vector3cd> H = {Eigen::Vector3cd::Zero()};
    std::vector<double> area = {1.0};
    std::vector<double> theta = {0.0, M_PI / 2};
    auto pat = stratton_chu_2d(r, n, E, H, area, theta, 0.0, 2 * M_PI);
    assert(pat.size() == 2);
    assert(std::abs(pat[0].e_theta) < 1e-12);
  }

  // Two-point surface with nonzero E and H; expect phi-polarized far field
  {
    std::vector<Eigen::Vector3d> r = {Eigen::Vector3d(0.5, 0.0, 0.0),
                                      Eigen::Vector3d(-0.5, 0.0, 0.0)};
    std::vector<Eigen::Vector3d> n(2, Eigen::Vector3d::UnitZ());
    std::vector<Eigen::Vector3cd> E(2, Eigen::Vector3cd::UnitX());
    std::vector<Eigen::Vector3cd> H(2, Eigen::Vector3cd::UnitY());
    std::vector<double> area = {1.0, 1.0};
    std::vector<double> theta = {0.0, M_PI / 2, M_PI};
    auto pat = stratton_chu_2d(r, n, E, H, area, theta, 0.0, 2 * M_PI);
    assert(pat.size() == 3);

    // Far field should be purely phi-polarized
    for (const auto &p : pat) {
      assert(std::abs(p.e_theta) < 1e-12);
    }

    double mag0 = std::abs(pat[0].e_phi);
    double mag1 = std::abs(pat[1].e_phi);
    double mag2 = std::abs(pat[2].e_phi);

    // Non-zero far field with symmetry about theta=pi/2
    assert(mag0 > 1e-3 && mag1 > 1e-3);
    assert(std::abs(mag0 - mag2) / mag0 < 1e-12);
    // Broadside maximum at theta=0 greater than at theta=pi/2
    assert(mag0 > 10 * mag1);
  }

  return 0;
}
