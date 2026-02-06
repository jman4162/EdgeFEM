#include "edgefem/array/embedded_pattern.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>

#include "edgefem/solver.hpp"

namespace edgefem {

namespace {

// Compute far-field point using simplified Huygens integration
// This is a placeholder - full implementation would use post/ntf.hpp
FFPoint2D compute_ff_point(const Eigen::VectorXcd &E_field, const Mesh &mesh,
                           double theta, double phi, double k0) {
  FFPoint2D point;
  point.theta = theta;

  // Simplified far-field computation
  // Real implementation would integrate tangential E and H over Huygens surface
  std::complex<double> E_ff(0, 0);

  // Approximate: sum contributions from near-field edges
  // weighted by their angular alignment with observation direction
  Eigen::Vector3d r_hat;
  r_hat.x() = std::sin(theta) * std::cos(phi);
  r_hat.y() = std::sin(theta) * std::sin(phi);
  r_hat.z() = std::cos(theta);

  for (size_t i = 0;
       i < mesh.edges.size() && i < static_cast<size_t>(E_field.size()); ++i) {
    const auto &edge = mesh.edges[i];
    Eigen::Vector3d p0 = mesh.nodes[mesh.nodeIndex.at(edge.n0)].xyz;
    Eigen::Vector3d p1 = mesh.nodes[mesh.nodeIndex.at(edge.n1)].xyz;
    Eigen::Vector3d edge_center = 0.5 * (p0 + p1);
    Eigen::Vector3d edge_vec = (p1 - p0).normalized();

    // Phase from observation point
    double phase_delay = k0 * edge_center.dot(r_hat);
    std::complex<double> phase_factor(std::cos(phase_delay),
                                      -std::sin(phase_delay));

    // Angular weighting (cross product with observation direction)
    double angular_weight = (edge_vec.cross(r_hat)).norm();

    E_ff += E_field(static_cast<int>(i)) * phase_factor * angular_weight;
  }

  point.magnitude = std::abs(E_ff);
  point.phase = std::arg(E_ff);

  return point;
}

} // namespace

std::vector<EmbeddedPattern>
compute_embedded_patterns(const Mesh &mesh, const MaxwellParams &params,
                          const BC &bc, const std::vector<WavePort> &ports,
                          const ArrayPatternConfig &config) {

  std::vector<EmbeddedPattern> patterns;
  patterns.reserve(ports.size());

  for (size_t i = 0; i < ports.size(); ++i) {
    patterns.push_back(compute_single_embedded_pattern(
        mesh, params, bc, ports, static_cast<int>(i), config));
  }

  return patterns;
}

EmbeddedPattern compute_single_embedded_pattern(
    const Mesh &mesh, const MaxwellParams &params, const BC &bc,
    const std::vector<WavePort> &ports, int active_port,
    const ArrayPatternConfig &config) {

  EmbeddedPattern pattern;
  pattern.port_index = active_port;
  pattern.frequency = params.omega / (2.0 * M_PI);

  // Estimate element position from port edge locations
  if (active_port >= 0 && active_port < static_cast<int>(ports.size())) {
    const auto &port = ports[active_port];
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int edge_idx : port.edges) {
      const auto &edge = mesh.edges[edge_idx];
      center += mesh.nodes[mesh.nodeIndex.at(edge.n0)].xyz;
      center += mesh.nodes[mesh.nodeIndex.at(edge.n1)].xyz;
    }
    if (!port.edges.empty()) {
      center /= (2.0 * port.edges.size());
    }
    pattern.element_position = center;
  }

  // Solve for E-field with active port excited
  auto asmbl = assemble_maxwell(mesh, params, bc, ports, active_port);
  auto result = solve_linear(asmbl.A, asmbl.b, {});

  if (!result.converged) {
    // Return empty pattern if solve failed
    return pattern;
  }

  // Compute E-plane pattern (phi = phi_e_plane)
  pattern.pattern_phi0.reserve(config.theta_rad.size());
  for (double theta : config.theta_rad) {
    pattern.pattern_phi0.push_back(
        compute_ff_point(result.x, mesh, theta, config.phi_e_plane, config.k0));
  }

  // Compute H-plane pattern (phi = phi_h_plane)
  pattern.pattern_phi90.reserve(config.theta_rad.size());
  for (double theta : config.theta_rad) {
    pattern.pattern_phi90.push_back(
        compute_ff_point(result.x, mesh, theta, config.phi_h_plane, config.k0));
  }

  // Normalize both patterns
  normalize_pattern(pattern.pattern_phi0);
  normalize_pattern(pattern.pattern_phi90);

  return pattern;
}

void normalize_pattern(std::vector<FFPoint2D> &pattern) {
  if (pattern.empty())
    return;

  // Find peak magnitude
  double peak = 0.0;
  for (const auto &pt : pattern) {
    if (pt.magnitude > peak) {
      peak = pt.magnitude;
    }
  }

  // Normalize to peak = 1.0
  if (peak > 1e-20) {
    for (auto &pt : pattern) {
      pt.magnitude /= peak;
    }
  }
}

} // namespace edgefem
