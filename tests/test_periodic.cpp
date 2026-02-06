#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include "edgefem/periodic.hpp"

using namespace edgefem;

void test_floquet_phase_from_angle_broadside() {
  std::cout << "test_floquet_phase_from_angle_broadside..." << std::endl;

  Eigen::Vector3d period(0.015, 0.0, 0.0);  // 15mm period in x
  double theta = 0.0;  // Broadside
  double phi = 0.0;
  double k0 = 2 * M_PI * 10e9 / 299792458.0;  // 10 GHz

  auto phase = floquet_phase_from_angle(period, theta, phi, k0);

  // At broadside, no phase shift
  assert(std::abs(phase - std::complex<double>(1.0, 0.0)) < 1e-10);

  std::cout << "  Broadside phase = " << phase << std::endl;
}

void test_floquet_phase_from_angle_scan() {
  std::cout << "test_floquet_phase_from_angle_scan..." << std::endl;

  Eigen::Vector3d period(0.015, 0.0, 0.0);  // 15mm period in x
  double theta = M_PI / 6.0;  // 30 degrees from broadside
  double phi = 0.0;  // E-plane
  double k0 = 2 * M_PI * 10e9 / 299792458.0;  // 10 GHz

  auto phase = floquet_phase_from_angle(period, theta, phi, k0);

  // kx = k0 * sin(30°) = k0 * 0.5
  // phase = exp(j * kx * Lx) = exp(j * k0 * 0.5 * 0.015)
  double kx = k0 * std::sin(theta);
  std::complex<double> expected(std::cos(kx * 0.015), std::sin(kx * 0.015));

  assert(std::abs(phase - expected) < 1e-10);

  std::cout << "  30° scan phase = " << phase << std::endl;
  std::cout << "  Phase angle = " << std::arg(phase) * 180 / M_PI << " deg"
            << std::endl;
}

void test_set_floquet_phase() {
  std::cout << "test_set_floquet_phase..." << std::endl;

  PeriodicBC pbc;
  pbc.period_vector = Eigen::Vector3d(0.015, 0.0, 0.0);

  double k0 = 2 * M_PI * 10e9 / 299792458.0;
  double theta = M_PI / 4.0;  // 45 degrees

  Eigen::Vector2d k_trans(k0 * std::sin(theta), 0.0);
  set_floquet_phase(pbc, k_trans);

  // Check phase was set correctly
  double expected_arg = k_trans.x() * pbc.period_vector.x();
  std::complex<double> expected(std::cos(expected_arg), std::sin(expected_arg));

  assert(std::abs(pbc.phase_shift - expected) < 1e-10);

  std::cout << "  Phase set to " << pbc.phase_shift << std::endl;
}

void test_periodic_pair_structure() {
  std::cout << "test_periodic_pair_structure..." << std::endl;

  PeriodicPair pair;
  pair.master_edge = 10;
  pair.slave_edge = 20;
  pair.master_orient = 1;
  pair.slave_orient = -1;
  pair.translation = Eigen::Vector3d(0.015, 0.0, 0.0);

  assert(pair.master_edge == 10);
  assert(pair.slave_edge == 20);
  assert(pair.master_orient * pair.slave_orient == -1);  // Opposite orientations

  std::cout << "  PeriodicPair structure valid" << std::endl;
}

void test_validate_periodic_bc_empty() {
  std::cout << "test_validate_periodic_bc_empty..." << std::endl;

  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");
  PeriodicBC pbc;  // Empty

  assert(!validate_periodic_bc(mesh, pbc));

  std::cout << "  Empty PBC correctly flagged as invalid" << std::endl;
}

void test_validate_periodic_bc_valid() {
  std::cout << "test_validate_periodic_bc_valid..." << std::endl;

  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");

  // Create a valid-looking PBC (may not correspond to actual mesh edges)
  PeriodicBC pbc;
  pbc.period_vector = Eigen::Vector3d(1.0, 0.0, 0.0);
  pbc.phase_shift = std::complex<double>(1.0, 0.0);

  // Add some valid pairs
  if (mesh.edges.size() >= 10) {
    PeriodicPair pair;
    pair.master_edge = 0;
    pair.slave_edge = 1;
    pair.master_orient = 1;
    pair.slave_orient = 1;
    pbc.pairs.push_back(pair);

    assert(validate_periodic_bc(mesh, pbc));
    std::cout << "  Valid PBC correctly validated" << std::endl;
  }
}

void test_count_surface_edges() {
  std::cout << "test_count_surface_edges..." << std::endl;

  Mesh mesh = load_gmsh_v2("examples/cube_cavity.msh");

  // Tag 1 is typically PEC boundary
  int count = count_surface_edges(mesh, 1);

  assert(count >= 0);
  std::cout << "  Surface with tag 1 has " << count << " edges" << std::endl;
}

void test_floquet_phase_grating_lobe() {
  std::cout << "test_floquet_phase_grating_lobe..." << std::endl;

  // At grating lobe onset: sin(theta) = lambda/d - 1
  // For d = lambda/2: sin(theta) = 2 - 1 = 1, theta = 90 deg

  double freq = 10e9;
  double c0 = 299792458.0;
  double lambda = c0 / freq;
  double period = lambda / 2.0;  // Half-wavelength spacing

  Eigen::Vector3d period_vec(period, 0.0, 0.0);
  double k0 = 2 * M_PI / lambda;
  double theta = M_PI / 2.0;  // 90 degrees (horizon)
  double phi = 0.0;

  auto phase = floquet_phase_from_angle(period_vec, theta, phi, k0);

  // At 90 deg: kx = k0, phase = exp(j * k0 * d) = exp(j * pi) = -1
  std::complex<double> expected(-1.0, 0.0);

  // Allow some numerical tolerance
  assert(std::abs(phase - expected) < 1e-6);

  std::cout << "  Phase at grating lobe onset = " << phase << std::endl;
}

int main() {
  test_floquet_phase_from_angle_broadside();
  test_floquet_phase_from_angle_scan();
  test_set_floquet_phase();
  test_periodic_pair_structure();
  test_validate_periodic_bc_empty();
  test_validate_periodic_bc_valid();
  test_count_surface_edges();
  test_floquet_phase_grating_lobe();
  std::cout << "All periodic BC tests passed!" << std::endl;
  return 0;
}
