#include "vectorem/io/touchstone.hpp"

#include <fstream>
#include <iomanip>

namespace vectorem {

void write_touchstone(const std::string &path, const std::vector<double> &freq,
                      const std::vector<SParams2> &data) {
  std::ofstream ofs(path);
  ofs << "# Hz S RI R 50\n";
  ofs << std::setprecision(12);
  for (size_t i = 0; i < freq.size(); ++i) {
    const auto &s = data[i];
    ofs << freq[i] << ' ' << s.s11.real() << ' ' << s.s11.imag() << ' '
        << s.s21.real() << ' ' << s.s21.imag() << ' ' << s.s12.real() << ' '
        << s.s12.imag() << ' ' << s.s22.real() << ' ' << s.s22.imag()
        << '\n';
  }
}

} // namespace vectorem

