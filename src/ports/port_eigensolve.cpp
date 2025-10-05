#include "vectorem/ports/port_eigensolve.hpp"

#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "vectorem/fem.hpp"

namespace vectorem {

namespace {
constexpr double c0 = 299792458.0;     // speed of light in vacuum (m/s)
constexpr double eta0 = 376.730313668; // impedance of free space (ohms)
constexpr double mu0 = 4.0 * M_PI * 1e-7;
constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

std::complex<double> ensure_sqrt(const std::complex<double> &value) {
  return std::sqrt(value);
}

Eigen::Matrix<double, 3, 2> shape_gradients(const Mesh &mesh,
                                            const Element &tri) {
  const auto &n0 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[0]));
  const auto &n1 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[1]));
  const auto &n2 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[2]));
  Eigen::Matrix3d C;
  C << 1.0, n0.xyz.x(), n0.xyz.y(), 1.0, n1.xyz.x(), n1.xyz.y(), 1.0,
      n2.xyz.x(), n2.xyz.y();
  return C.inverse().block<3, 2>(0, 1);
}

double tri_area(const Mesh &mesh, const Element &tri) {
  const auto &n0 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[0]));
  const auto &n1 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[1]));
  const auto &n2 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[2]));
  Eigen::Matrix3d C;
  C << 1.0, n0.xyz.x(), n0.xyz.y(), 1.0, n1.xyz.x(), n1.xyz.y(), 1.0,
      n2.xyz.x(), n2.xyz.y();
  return 0.5 * std::abs(C.determinant());
}
} // namespace

PortMode solve_te10_mode(const RectWaveguidePort &port, double freq) {
  PortMode mode{};
  mode.pol = ModePolarization::TE;
  mode.fc = c0 / (2.0 * port.a);
  const double k = 2.0 * M_PI * freq / c0;
  const double kc = M_PI / port.a; // TE10 cutoff wavenumber
  mode.kc = kc;
  mode.omega = 2.0 * M_PI * freq;
  mode.mu = mu0;
  mode.eps = eps0;
  const double beta_real = std::sqrt(std::max(0.0, k * k - kc * kc));
  mode.beta = beta_real;
  if (beta_real > 0) {
    mode.Z0 = eta0 * k / beta_real;
  } else {
    mode.Z0 = 0.0;
  }
  return mode;
}

SParams2 straight_waveguide_sparams(const RectWaveguidePort &port,
                                    double length, double freq) {
  const double k = 2.0 * M_PI * freq / c0;
  const double kc = M_PI / port.a;
  SParams2 s{};
  if (k <= kc) {
    // Below cutoff: mode is evanescent, full reflection
    s.s11 = 1.0;
    s.s22 = 1.0;
    s.s21 = 0.0;
    s.s12 = 0.0;
  } else {
    const double beta = std::sqrt(k * k - kc * kc);
    const std::complex<double> phase =
        std::exp(std::complex<double>(0.0, -beta * length));
    s.s11 = 0.0;
    s.s22 = 0.0;
    s.s21 = phase;
    s.s12 = phase;
  }
  return s;
}

} // namespace vectorem

namespace vectorem {

using TripletC = Eigen::Triplet<std::complex<double>>;

std::vector<PortMode>
solve_port_eigens(const Mesh &mesh, int num_modes, double omega,
                  std::complex<double> eps_r, std::complex<double> mu_r,
                  ModePolarization pol) {
  if (mesh.tris.empty()) {
    throw std::runtime_error("Port eigensolver requires a 2D mesh.");
  }

  if (omega == 0.0) {
    throw std::runtime_error("Port eigensolver requires non-zero frequency.");
  }

  const int num_nodes = mesh.nodes.size();

  // Identify PEC nodes and create a map for non-PEC ("free") DOFs.
  std::vector<bool> is_pec_node(num_nodes, false);
  if (pol == ModePolarization::TM) {
    const int pec_tag = 1;
    for (const auto &b : mesh.boundary_lines) {
      if (b.phys == pec_tag) {
        is_pec_node.at(mesh.nodeIndex.at(b.n0)) = true;
        is_pec_node.at(mesh.nodeIndex.at(b.n1)) = true;
      }
    }
  }

  std::vector<int> node_to_dof(num_nodes, -1);
  int num_dofs = 0;
  for (int i = 0; i < num_nodes; ++i) {
    if (!is_pec_node[i]) {
      node_to_dof[i] = num_dofs++;
    }
  }

  std::vector<TripletC> A_triplets, B_triplets;
  A_triplets.reserve(mesh.tris.size() * 9);
  B_triplets.reserve(mesh.tris.size() * 9);

  for (const auto &tri : mesh.tris) {
    const auto &n0 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[0]));
    const auto &n1 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[1]));
    const auto &n2 = mesh.nodes.at(mesh.nodeIndex.at(tri.conn[2]));

    Eigen::Matrix3d C;
    C << 1, n0.xyz.x(), n0.xyz.y(), 1, n1.xyz.x(), n1.xyz.y(), 1, n2.xyz.x(), n2.xyz.y();

    double area = 0.5 * std::abs(C.determinant());

    Eigen::Matrix<double, 3, 2> G = C.inverse().block<3, 2>(0, 1);
    Eigen::Matrix3d ke = area * G * G.transpose();
    Eigen::Matrix3d me = (area / 12.0) * Eigen::Matrix3d{{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};

    for (int i = 0; i < 3; ++i) {
      int node_i = mesh.nodeIndex.at(tri.conn[i]);
      int dof_i = node_to_dof[node_i];
      if (dof_i < 0) continue; // Skip PEC node

      for (int j = 0; j < 3; ++j) {
        int node_j = mesh.nodeIndex.at(tri.conn[j]);
        int dof_j = node_to_dof[node_j];
        if (dof_j < 0) continue; // Skip PEC node

        A_triplets.emplace_back(dof_i, dof_j, ke(i, j));
        B_triplets.emplace_back(dof_i, dof_j, me(i, j));
      }
    }
  }

  SpMatC A(num_dofs, num_dofs);
  SpMatC B(num_dofs, num_dofs);
  A.setFromTriplets(A_triplets.begin(), A_triplets.end());
  B.setFromTriplets(B_triplets.begin(), B_triplets.end());

  Eigen::MatrixXcd Ad = A;
  Eigen::MatrixXcd Bd = B;

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> ges(Ad, Bd);
  if (ges.info() != Eigen::Success) {
    throw std::runtime_error("Eigenvalue computation failed.");
  }

  std::vector<PortMode> modes;
  const auto &eigenvalues = ges.eigenvalues();
  const auto &eigenvectors = ges.eigenvectors();

  const std::complex<double> eps = eps0 * eps_r;
  const std::complex<double> mu = mu0 * mu_r;
  const std::complex<double> k = omega * ensure_sqrt(mu * eps);

  for (int i = 0; i < eigenvalues.size(); ++i) {
    double kc2 = eigenvalues[i];
    if (!std::isfinite(kc2) || kc2 < 1e-12)
      continue;
    double kc = std::sqrt(kc2);
    double fc = kc * c0 / (2.0 * M_PI);

    PortMode mode;
    mode.pol = pol;
    mode.fc = fc;
    mode.kc = kc;
    mode.omega = omega;
    mode.eps = eps;
    mode.mu = mu;

    std::complex<double> beta = ensure_sqrt(k * k - kc2);
    mode.beta = beta;

    if (pol == ModePolarization::TE) {
      mode.Z0 = (beta == std::complex<double>(0.0))
                    ? std::complex<double>(0.0)
                    : (omega * mu / beta);
    } else {
      mode.Z0 = (beta == std::complex<double>(0.0))
                    ? std::complex<double>(0.0)
                    : (beta / (omega * eps));
    }

    Eigen::VectorXcd vec = eigenvectors.col(i);
    Eigen::VectorXcd field_full = Eigen::VectorXcd::Zero(num_nodes);
    for (int n = 0; n < num_nodes; ++n) {
      int dof = node_to_dof[n];
      if (dof >= 0)
        field_full[n] = vec[dof];
    }

    std::complex<double> power = 0.0;
    const std::complex<double> j(0.0, 1.0);
    for (const auto &tri : mesh.tris) {
      Eigen::Matrix<double, 3, 2> G = shape_gradients(mesh, tri);
      Eigen::Vector3cd values;
      for (int k = 0; k < 3; ++k) {
        values[k] =
            field_full[mesh.nodeIndex.at(tri.conn[k])];
      }
      Eigen::Vector2cd grad = G.transpose() * values;

      Eigen::Matrix<std::complex<double>, 3, 1> Et;
      Eigen::Matrix<std::complex<double>, 3, 1> Ht;
      if (pol == ModePolarization::TE) {
        auto factor_e = j * omega * mu / (kc * kc);
        auto factor_h = beta / (kc * kc);
        Et << -factor_e * grad(1), factor_e * grad(0), 0.0;
        Ht << factor_h * grad(0), factor_h * grad(1), 0.0;
      } else {
        auto factor_e = -beta / (kc * kc);
        auto factor_h = 1.0 / (j * omega * mu);
        Et << factor_e * grad(0), factor_e * grad(1), 0.0;
        Ht << -factor_h * grad(1), factor_h * grad(0), 0.0;
      }

      std::complex<double> p =
          Et(0) * std::conj(Ht(1)) - Et(1) * std::conj(Ht(0));
      power += 0.5 * tri_area(mesh, tri) * p;
    }

    double power_real = std::real(power);
    if (power_real <= 0.0)
      power_real = std::abs(power);
    if (power_real <= 0.0)
      continue;

    double scale = 1.0 / std::sqrt(power_real);
    field_full *= scale;
    mode.field = field_full;
    modes.push_back(std::move(mode));
    if (static_cast<int>(modes.size()) >= num_modes)
      break;
  }

  std::sort(modes.begin(), modes.end(),
            [](const PortMode &a, const PortMode &b) { return a.fc < b.fc; });

  if (static_cast<int>(modes.size()) > num_modes)
    modes.resize(num_modes);

  return modes;
}

}
