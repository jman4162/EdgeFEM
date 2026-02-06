#pragma once

#include <complex>
#include <vector>

#include <Eigen/Core>

#include "edgefem/mesh.hpp"

namespace edgefem {

/// Pair of edges linked by periodic boundary condition.
struct PeriodicPair {
  int master_edge;     // Global edge index on master surface
  int slave_edge;      // Global edge index on slave surface
  int master_orient;   // Orientation of edge on master (+1 or -1)
  int slave_orient;    // Orientation of edge on slave (+1 or -1)
  Eigen::Vector3d translation;  // Vector from master to slave
};

/// Periodic boundary condition specification.
struct PeriodicBC {
  std::vector<PeriodicPair> pairs;  // Master-slave edge pairs
  Eigen::Vector3d period_vector;    // Lattice translation vector
  std::complex<double> phase_shift; // Floquet phase: exp(j*k·L)
};

/// Build periodic edge pairs from mesh surface tags.
/// Matches edges on the slave surface to corresponding edges on the master
/// surface based on their geometric positions translated by the period vector.
///
/// @param mesh Volume mesh containing both surfaces
/// @param master_tag Physical tag of master surface
/// @param slave_tag Physical tag of slave surface
/// @param period_vector Translation from master to slave surface
/// @param tolerance Position matching tolerance
/// @return PeriodicBC structure with edge pairs and period vector
PeriodicBC build_periodic_pairs(const Mesh &mesh, int master_tag, int slave_tag,
                                const Eigen::Vector3d &period_vector,
                                double tolerance = 1e-9);

/// Validate that periodic BC pairs are properly matched.
/// Checks that all edges on both surfaces are paired and orientations are
/// consistent.
///
/// @param mesh Volume mesh
/// @param pbc Periodic BC to validate
/// @return true if valid, false if issues found
bool validate_periodic_bc(const Mesh &mesh, const PeriodicBC &pbc);

/// Set Floquet phase shift for given transverse wave vector.
/// The phase shift is computed as exp(j * k_transverse · period_vector).
///
/// @param pbc PeriodicBC to modify
/// @param k_transverse Transverse wave vector [kx, ky] (rad/m)
void set_floquet_phase(PeriodicBC &pbc, const Eigen::Vector2d &k_transverse);

/// Compute Floquet phase shift from incidence angles.
/// @param period_vector Lattice translation vector
/// @param theta Elevation angle from z-axis (radians)
/// @param phi Azimuth angle from x-axis (radians)
/// @param k0 Free-space wavenumber (rad/m)
/// @return Phase shift exp(j * k_inc · period_vector)
std::complex<double> floquet_phase_from_angle(
    const Eigen::Vector3d &period_vector, double theta, double phi, double k0);

/// Count edges on a surface with given physical tag.
/// Useful for verifying periodic BC matching.
///
/// @param mesh Volume mesh
/// @param surface_tag Physical tag of surface
/// @return Number of edges on the surface
int count_surface_edges(const Mesh &mesh, int surface_tag);

} // namespace edgefem
