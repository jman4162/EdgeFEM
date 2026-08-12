#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// assert() vanishes under -DNDEBUG (release builds); the dipole test below
// must fail loudly regardless of build type.
#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::fprintf(stderr, "CHECK failed: %s (%s:%d)\n", #cond, __FILE__,      \
                   __LINE__);                                                  \
      std::exit(1);                                                            \
    }                                                                          \
  } while (0)

#include "edgefem/post/ntf.hpp"

using namespace edgefem;

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

    // Hertzian dipole absolute test: exact near fields sampled on a Huygens
  // sphere must reproduce the analytic far field in shape, polarization,
  // AND absolute magnitude. This is the test that arbitrates the kernel:
  // the pre-fix kernel (jk0 on the M term only) fails it by orders of
  // magnitude and with a distorted pattern.
  //
  // z-directed dipole, moment I*dl = 1 A.m, e^{-jwt} convention (outgoing
  // e^{+jkr}), exact fields (Balanis 4-8a..4-10c with conjugated j):
  //   E_r     = (Z0 Idl cos(th) / (2 pi r^2)) (1 + 1/(jkr) ... ) etc.
  // Far field: |E_theta| = Z0 k Idl sin(th) / (4 pi r)  -> r-normalized
  // amplitude Z0 k Idl sin(th) / (4 pi).
  {
    const double k0 = 2.0 * M_PI;         // 1 m wavelength
    const double a = 0.75;                // Huygens sphere radius (ka ~ 4.7)
    const double Z0l = 376.730313668;
    const std::complex<double> j(0.0, 1.0);

    // Sample sphere: theta x phi quadrature (midpoint in both angles)
    const int Ns_th = 90, Ns_ph = 180;
    std::vector<Eigen::Vector3d> r, n;
    std::vector<Eigen::Vector3cd> E, H;
    std::vector<double> area;
    for (int it = 0; it < Ns_th; ++it) {
      double th = (it + 0.5) * M_PI / Ns_th;
      for (int ip = 0; ip < Ns_ph; ++ip) {
        double ph = (ip + 0.5) * 2.0 * M_PI / Ns_ph;
        Eigen::Vector3d rhat(std::sin(th) * std::cos(ph),
                             std::sin(th) * std::sin(ph), std::cos(th));
        Eigen::Vector3d th_hat(std::cos(th) * std::cos(ph),
                               std::cos(th) * std::sin(ph), -std::sin(th));
        Eigen::Vector3d ph_hat(-std::sin(ph), std::cos(ph), 0.0);

        // Exact dipole fields at radius a, written in Balanis's e^{+jwt}
        // convention (4-8a..4-10c, outgoing e^{-jka}, u = 1/(jka)), then
        // conjugated wholesale to the kernel's e^{-jwt} convention —
        // conjugating term-by-term is where sign errors breed.
        std::complex<double> u = 1.0 / (j * k0 * a);
        std::complex<double> emjka = std::exp(-j * k0 * a);
        std::complex<double> Er =
            Z0l * std::cos(th) / (2.0 * M_PI * a * a) * (1.0 + u) * emjka;
        std::complex<double> Eth = j * Z0l * k0 * std::sin(th) /
                                   (4.0 * M_PI * a) * (1.0 + u + u * u) * emjka;
        std::complex<double> Hph =
            j * k0 * std::sin(th) / (4.0 * M_PI * a) * (1.0 + u) * emjka;

        r.push_back(a * rhat);
        n.push_back(rhat);
        Eigen::Vector3cd Evec = Er * rhat.cast<std::complex<double>>() +
                                Eth * th_hat.cast<std::complex<double>>();
        Eigen::Vector3cd Hvec = Hph * ph_hat.cast<std::complex<double>>();
        E.push_back(Evec.conjugate());
        H.push_back(Hvec.conjugate());
        double dA = a * a * std::sin(th) * (M_PI / Ns_th) * (2.0 * M_PI / Ns_ph);
        area.push_back(dA);
      }
    }

    // Far-field grid
    const int Nth = 19, Nph = 4;
    std::vector<double> theta_g(Nth), phi_g(Nph);
    for (int i = 0; i < Nth; ++i)
      theta_g[i] = i * M_PI / (Nth - 1);
    for (int i = 0; i < Nph; ++i)
      phi_g[i] = i * M_PI / 2.0;

    FFPattern3D pat = stratton_chu_3d(r, n, E, H, area, theta_g, phi_g, k0);

    // Expected r-normalized far-field magnitude: Z0 k0 sin(th) / (4 pi)
    // (I dl = 1). Check absolute magnitude at several angles.
    for (int it = 1; it + 1 < Nth; ++it) { // skip poles (sin=0)
      double th = theta_g[it];
      double expected = Z0l * k0 * std::sin(th) / (4.0 * M_PI);
      for (int ip = 0; ip < Nph; ++ip) {
        double got = std::abs(pat.E_theta(it, ip));
        CHECK(std::abs(got - expected) / expected < 0.02); // 2% quadrature
        CHECK(std::abs(pat.E_phi(it, ip)) < 1e-3 * expected); // pure theta-pol
      }
    }

    // Directivity of a Hertzian dipole is exactly 1.5 (1.76 dBi)
    {
      std::vector<double> theta_d(61), phi_d(121);
      for (int i = 0; i < 61; ++i) theta_d[i] = i * M_PI / 60;
      for (int i = 0; i < 121; ++i) phi_d[i] = i * 2.0 * M_PI / 120;
      FFPattern3D pat_d = stratton_chu_3d(r, n, E, H, area, theta_d, phi_d, k0);
      double D = compute_directivity(pat_d);
      CHECK(std::abs(D - 1.5) < 0.02);
    }
  }

  return 0;
}
