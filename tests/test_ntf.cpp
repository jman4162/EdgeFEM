#include "test_check.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

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
    CHECK(pat.size() == 2);
    CHECK(std::abs(pat[0].e_theta) < 1e-12);
  }

  // stratton_chu_2d is the phi = const cut of stratton_chu_3d.
  {
    std::vector<Eigen::Vector3d> r = {Eigen::Vector3d(0.3, -0.1, 0.0),
                                      Eigen::Vector3d(-0.2, 0.4, 0.1)};
    std::vector<Eigen::Vector3d> n = {Eigen::Vector3d::UnitZ(),
                                      Eigen::Vector3d::UnitX()};
    std::vector<Eigen::Vector3cd> E = {
        Eigen::Vector3cd(1.0, std::complex<double>(0.0, 0.5), 0.0),
        Eigen::Vector3cd(0.0, 0.2, std::complex<double>(-0.3, 0.1))};
    std::vector<Eigen::Vector3cd> H = {
        Eigen::Vector3cd(0.0, 1e-3, 0.0),
        Eigen::Vector3cd(0.0, std::complex<double>(0.0, 2e-3), 1e-3)};
    std::vector<double> area = {0.5, 0.25};
    std::vector<double> theta = {0.0, 0.4, M_PI / 2, 2.0, M_PI};
    const double phi = 0.7, k0 = 2 * M_PI;
    auto cut = stratton_chu_2d(r, n, E, H, area, theta, phi, k0);
    FFPattern3D full = stratton_chu_3d(r, n, E, H, area, theta, {phi}, k0);
    CHECK(cut.size() == theta.size());
    for (size_t i = 0; i < theta.size(); ++i) {
      CHECK(std::abs(cut[i].e_theta - full.E_theta(i, 0)) < 1e-12);
      CHECK(std::abs(cut[i].e_phi - full.E_phi(i, 0)) < 1e-12);
    }
  }

  // exp(+j omega t) aperture: a plane wave leaving z = 0 toward
  // k = k0 (sin th0, 0, cos th0), E = x exp(-j k0 sin(th0) x),
  // H = (k x E) / Z0, on a 6 x 6 wavelength patch with outward n = +z.
  // Its far field peaks at (theta, phi) = (th0, 0). Behind the aperture, at
  // (pi - th0, 0), the array factor is the same and the Huygens element
  // factor 1 + cos(th0) cos(theta) gives the power ratio
  // ((1 - cos^2 th0) / (1 + cos^2 th0))^2 = 0.0204 for th0 = 30 deg. A
  // kernel with the exp(-j omega t) phase puts the peak at phi = pi.
  {
    const double k0 = 2.0 * M_PI, Z0l = 376.730313668;
    const double th0 = M_PI / 6.0;
    const std::complex<double> j(0.0, 1.0);
    const int Na = 60;
    const double side = 6.0, h = side / Na;
    std::vector<Eigen::Vector3d> r, n;
    std::vector<Eigen::Vector3cd> E, H;
    std::vector<double> area;
    for (int ix = 0; ix < Na; ++ix) {
      for (int iy = 0; iy < Na; ++iy) {
        const double x = -side / 2 + (ix + 0.5) * h;
        const double y = -side / 2 + (iy + 0.5) * h;
        const std::complex<double> ph = std::exp(-j * k0 * std::sin(th0) * x);
        r.emplace_back(x, y, 0.0);
        n.push_back(Eigen::Vector3d::UnitZ());
        E.push_back(Eigen::Vector3cd(ph, 0.0, 0.0));
        H.push_back(Eigen::Vector3cd(0.0, std::cos(th0) * ph / Z0l, 0.0));
        area.push_back(h * h);
      }
    }
    std::vector<double> theta;
    for (int i = 0; i <= 180; ++i)
      theta.push_back(i * M_PI / 180.0);
    FFPattern3D pat = stratton_chu_3d(r, n, E, H, area, theta, {0.0, M_PI}, k0);
    Eigen::MatrixXd pwr = pat.power_pattern();
    Eigen::Index i_max = 0, j_max = 0;
    pwr.maxCoeff(&i_max, &j_max);
    CHECK(j_max == 0);                                        // phi = 0 side
    CHECK(std::abs(theta[i_max] - th0) < 2.0 * M_PI / 180.0); // at th0
    const double peak = pwr(i_max, j_max);
    CHECK(pwr(30, 1) < 1e-2 * peak); // mirror direction (th0, phi = pi)
    const double c2 = std::cos(th0) * std::cos(th0);
    const double back = std::pow((1.0 - c2) / (1.0 + c2), 2);
    CHECK(std::abs(pwr(150, 0) / peak - back) < 0.1 * back);
  }

  // Hertzian dipole absolute test: exact near fields sampled on a Huygens
  // sphere must reproduce the analytic far field in shape, polarization,
  // AND absolute magnitude. This is the test that arbitrates the kernel:
  // the pre-fix kernel (jk0 on the M term only) fails it by orders of
  // magnitude and with a distorted pattern.
  //
  // z-directed dipole, moment I*dl = 1 A.m, exp(+j omega t) convention
  // (outgoing exp(-jkr)), exact fields from Balanis 4-8a..4-10c. Far field:
  // E_theta = j Z0 k Idl sin(th) / (4 pi r) exp(-jkr); the r-normalized
  // pattern (exp(-jkr)/r removed) is j Z0 k sin(th) / (4 pi).
  {
    const double k0 = 2.0 * M_PI; // 1 m wavelength
    const double a = 0.75;        // Huygens sphere radius (ka ~ 4.7)
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

        // Exact dipole fields at radius a (Balanis 4-8a..4-10c, outgoing
        // exp(-jka), u = 1/(jka)), in the solver's exp(+j omega t)
        // convention.
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
        E.push_back(Evec);
        H.push_back(Hvec);
        double dA =
            a * a * std::sin(th) * (M_PI / Ns_th) * (2.0 * M_PI / Ns_ph);
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

    // Expected r-normalized far field: j Z0 k0 sin(th) / (4 pi) (I dl = 1),
    // magnitude and phase, at several angles.
    for (int it = 1; it + 1 < Nth; ++it) { // skip poles (sin=0)
      double th = theta_g[it];
      const std::complex<double> expected =
          j * Z0l * k0 * std::sin(th) / (4.0 * M_PI);
      for (int ip = 0; ip < Nph; ++ip) {
        const std::complex<double> got = pat.E_theta(it, ip);
        CHECK(std::abs(got - expected) / std::abs(expected) < 0.02);
        CHECK(std::abs(pat.E_phi(it, ip)) < 1e-3 * std::abs(expected));
      }
    }

    // Directivity of a Hertzian dipole is exactly 1.5 (1.76 dBi)
    {
      std::vector<double> theta_d(61), phi_d(121);
      for (int i = 0; i < 61; ++i)
        theta_d[i] = i * M_PI / 60;
      for (int i = 0; i < 121; ++i)
        phi_d[i] = i * 2.0 * M_PI / 120;
      FFPattern3D pat_d = stratton_chu_3d(r, n, E, H, area, theta_d, phi_d, k0);
      double D = compute_directivity(pat_d);
      CHECK(std::abs(D - 1.5) < 0.02);
    }
  }

  return edgefem_test::exit_code();
}
