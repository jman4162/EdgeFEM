#pragma once

#include <string>
#include <vector>

#include "vectorem/ports/port_eigensolve.hpp"

namespace vectorem {

void write_touchstone(const std::string &path, const std::vector<double> &freq,
                      const std::vector<SParams2> &data);

} // namespace vectorem
