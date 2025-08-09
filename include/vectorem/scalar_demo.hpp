#pragma once

#include <Eigen/Core>
#include <iostream>

namespace vectorem {

inline void run_scalar_demo() {
    Eigen::Vector3d v = Eigen::Vector3d::Ones();
    std::cout << "VectorEM scalar demo: " << v.transpose() << std::endl;
}

} // namespace vectorem

