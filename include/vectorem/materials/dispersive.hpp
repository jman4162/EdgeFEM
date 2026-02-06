#pragma once

#include <complex>
#include <memory>
#include <stdexcept>
#include <vector>

namespace vectorem {
namespace materials {

/// Base class for frequency-dependent (dispersive) material models.
/// Subclasses implement eval_eps() and optionally eval_mu() to return
/// the complex permittivity/permeability at a given angular frequency.
class DispersiveMaterial {
public:
  virtual ~DispersiveMaterial() = default;

  /// Evaluate complex relative permittivity at angular frequency omega (rad/s).
  /// @param omega Angular frequency (rad/s)
  /// @return Complex relative permittivity eps_r(omega)
  virtual std::complex<double> eval_eps(double omega) const = 0;

  /// Evaluate complex relative permeability at angular frequency omega (rad/s).
  /// Default implementation returns 1.0 (non-magnetic material).
  /// @param omega Angular frequency (rad/s)
  /// @return Complex relative permeability mu_r(omega)
  virtual std::complex<double> eval_mu(double omega) const {
    (void)omega;
    return std::complex<double>(1.0, 0.0);
  }
};

/// Debye relaxation model for polar dielectrics.
///
/// The Debye model describes frequency-dependent polarization with a single
/// relaxation time constant. It is commonly used for polar liquids (water),
/// biological tissues, and polymers.
///
/// Formula: eps(omega) = eps_inf + (eps_s - eps_inf) / (1 + j*omega*tau)
///
/// where:
///   eps_s   = static (DC) permittivity
///   eps_inf = high-frequency (optical) permittivity
///   tau     = relaxation time constant (seconds)
///
/// Example (water at 20°C):
///   eps_s = 80.1, eps_inf = 4.9, tau = 9.3e-12 s
class DebyeMaterial : public DispersiveMaterial {
public:
  /// Construct a Debye material model.
  /// @param eps_static Static (DC) relative permittivity
  /// @param eps_inf High-frequency relative permittivity
  /// @param tau Relaxation time constant in seconds (must be > 0)
  /// @throws std::invalid_argument if tau <= 0
  DebyeMaterial(double eps_static, double eps_inf, double tau)
      : eps_s_(eps_static), eps_inf_(eps_inf), tau_(tau) {
    if (tau <= 0.0) {
      throw std::invalid_argument("DebyeMaterial: tau must be positive");
    }
  }

  std::complex<double> eval_eps(double omega) const override {
    // eps(omega) = eps_inf + (eps_s - eps_inf) / (1 + j*omega*tau)
    std::complex<double> denom(1.0, omega * tau_);
    return eps_inf_ + (eps_s_ - eps_inf_) / denom;
  }

  double eps_static() const { return eps_s_; }
  double eps_inf() const { return eps_inf_; }
  double tau() const { return tau_; }

private:
  double eps_s_;    ///< Static (DC) permittivity
  double eps_inf_;  ///< High-frequency permittivity
  double tau_;      ///< Relaxation time (seconds)
};

/// Multi-pole Lorentz oscillator model for resonant materials.
///
/// The Lorentz model describes materials with atomic/molecular resonances,
/// such as dielectrics near absorption bands, metamaterials, and crystals.
///
/// Formula: eps(omega) = eps_inf + Sum_k [ delta_eps_k * omega0_k^2 /
///                                         (omega0_k^2 - omega^2 + j*gamma_k*omega) ]
///
/// Each pole is characterized by:
///   delta_eps = oscillator strength (contribution to permittivity)
///   omega0    = resonance angular frequency (rad/s)
///   gamma     = damping/loss factor (rad/s)
///
/// Example (single-pole resonance):
///   eps_inf = 1.0, delta_eps = 5.0, omega0 = 2*pi*10e9, gamma = 2*pi*0.5e9
class LorentzMaterial : public DispersiveMaterial {
public:
  /// Pole parameters for a single Lorentz oscillator.
  struct Pole {
    double delta_eps;  ///< Oscillator strength
    double omega0;     ///< Resonance angular frequency (rad/s)
    double gamma;      ///< Damping factor (rad/s)
  };

  /// Construct an empty Lorentz material with eps_inf = 1.0.
  /// Add poles using add_pole().
  LorentzMaterial() : eps_inf_(1.0) {}

  /// Construct a Lorentz material with specified high-frequency permittivity.
  /// @param eps_inf High-frequency permittivity
  explicit LorentzMaterial(double eps_inf) : eps_inf_(eps_inf) {}

  /// Add a Lorentz pole to the model.
  /// @param delta_eps Oscillator strength (permittivity contribution)
  /// @param omega0 Resonance angular frequency (rad/s, must be > 0)
  /// @param gamma Damping factor (rad/s, must be >= 0)
  /// @throws std::invalid_argument if omega0 <= 0 or gamma < 0
  void add_pole(double delta_eps, double omega0, double gamma) {
    if (omega0 <= 0.0) {
      throw std::invalid_argument("LorentzMaterial: omega0 must be positive");
    }
    if (gamma < 0.0) {
      throw std::invalid_argument("LorentzMaterial: gamma must be non-negative");
    }
    poles_.push_back({delta_eps, omega0, gamma});
  }

  std::complex<double> eval_eps(double omega) const override {
    std::complex<double> eps = eps_inf_;
    const double omega_sq = omega * omega;

    for (const auto &pole : poles_) {
      // delta_eps * omega0^2 / (omega0^2 - omega^2 + j*gamma*omega)
      double omega0_sq = pole.omega0 * pole.omega0;
      std::complex<double> denom(omega0_sq - omega_sq, pole.gamma * omega);
      eps += pole.delta_eps * omega0_sq / denom;
    }
    return eps;
  }

  double eps_inf() const { return eps_inf_; }
  void set_eps_inf(double eps_inf) { eps_inf_ = eps_inf; }
  const std::vector<Pole> &poles() const { return poles_; }
  size_t num_poles() const { return poles_.size(); }

private:
  double eps_inf_;          ///< High-frequency permittivity
  std::vector<Pole> poles_; ///< List of Lorentz poles
};

/// Drude free-electron model for metals and plasmas.
///
/// The Drude model describes the optical response of free electrons in metals
/// and plasmas. It produces negative permittivity below the plasma frequency.
///
/// Formula: eps(omega) = 1 - omega_p^2 / (omega^2 + j*gamma*omega)
///
/// where:
///   omega_p = plasma frequency (rad/s)
///   gamma   = collision/damping frequency (rad/s)
///
/// Example (gold at optical frequencies):
///   omega_p = 1.37e16 rad/s, gamma = 4.05e13 rad/s
///
/// Note: Real metals often require a Drude-Lorentz model combining both
/// free-electron (Drude) and bound-electron (Lorentz) contributions.
class DrudeMaterial : public DispersiveMaterial {
public:
  /// Construct a Drude material model.
  /// @param omega_p Plasma frequency in rad/s (must be > 0)
  /// @param gamma Collision/damping frequency in rad/s (must be >= 0)
  /// @throws std::invalid_argument if omega_p <= 0 or gamma < 0
  DrudeMaterial(double omega_p, double gamma)
      : omega_p_(omega_p), gamma_(gamma) {
    if (omega_p <= 0.0) {
      throw std::invalid_argument("DrudeMaterial: omega_p must be positive");
    }
    if (gamma < 0.0) {
      throw std::invalid_argument("DrudeMaterial: gamma must be non-negative");
    }
  }

  std::complex<double> eval_eps(double omega) const override {
    // eps(omega) = 1 - omega_p^2 / (omega^2 + j*gamma*omega)
    //            = 1 - omega_p^2 / (omega * (omega + j*gamma))
    if (omega == 0.0) {
      // DC limit: eps -> -infinity for pure Drude
      // Return a large negative value to indicate metallic behavior
      return std::complex<double>(-1e30, 0.0);
    }
    double omega_sq = omega * omega;
    double omega_p_sq = omega_p_ * omega_p_;
    std::complex<double> denom(omega_sq, gamma_ * omega);
    return std::complex<double>(1.0, 0.0) - omega_p_sq / denom;
  }

  double omega_p() const { return omega_p_; }
  double gamma() const { return gamma_; }

private:
  double omega_p_;  ///< Plasma frequency (rad/s)
  double gamma_;    ///< Collision/damping frequency (rad/s)
};

/// Combined Drude-Lorentz model for realistic metals.
///
/// This model combines free-electron (Drude) and bound-electron (Lorentz)
/// contributions, which is necessary for accurate modeling of real metals
/// across a wide frequency range.
///
/// Formula: eps(omega) = eps_inf - omega_p^2/(omega^2 + j*gamma_d*omega)
///                     + Sum_k [ delta_eps_k * omega0_k^2 /
///                               (omega0_k^2 - omega^2 + j*gamma_k*omega) ]
class DrudeLorentzMaterial : public DispersiveMaterial {
public:
  /// Lorentz pole parameters
  struct LorentzPole {
    double delta_eps;
    double omega0;
    double gamma;
  };

  /// Construct a Drude-Lorentz material.
  /// @param eps_inf High-frequency permittivity (background)
  /// @param omega_p Plasma frequency (rad/s)
  /// @param gamma_d Drude damping (rad/s)
  DrudeLorentzMaterial(double eps_inf, double omega_p, double gamma_d)
      : eps_inf_(eps_inf), omega_p_(omega_p), gamma_d_(gamma_d) {
    if (omega_p <= 0.0) {
      throw std::invalid_argument("DrudeLorentzMaterial: omega_p must be positive");
    }
    if (gamma_d < 0.0) {
      throw std::invalid_argument("DrudeLorentzMaterial: gamma_d must be non-negative");
    }
  }

  /// Add a Lorentz pole for interband transitions.
  void add_lorentz_pole(double delta_eps, double omega0, double gamma) {
    if (omega0 <= 0.0) {
      throw std::invalid_argument("DrudeLorentzMaterial: omega0 must be positive");
    }
    if (gamma < 0.0) {
      throw std::invalid_argument("DrudeLorentzMaterial: gamma must be non-negative");
    }
    lorentz_poles_.push_back({delta_eps, omega0, gamma});
  }

  std::complex<double> eval_eps(double omega) const override {
    std::complex<double> eps = eps_inf_;

    // Drude contribution
    if (omega != 0.0) {
      double omega_sq = omega * omega;
      double omega_p_sq = omega_p_ * omega_p_;
      std::complex<double> drude_denom(omega_sq, gamma_d_ * omega);
      eps -= omega_p_sq / drude_denom;
    } else {
      eps = std::complex<double>(-1e30, 0.0);
    }

    // Lorentz contributions
    double omega_sq = omega * omega;
    for (const auto &pole : lorentz_poles_) {
      double omega0_sq = pole.omega0 * pole.omega0;
      std::complex<double> denom(omega0_sq - omega_sq, pole.gamma * omega);
      eps += pole.delta_eps * omega0_sq / denom;
    }

    return eps;
  }

  double eps_inf() const { return eps_inf_; }
  double omega_p() const { return omega_p_; }
  double gamma_d() const { return gamma_d_; }
  const std::vector<LorentzPole> &lorentz_poles() const { return lorentz_poles_; }

private:
  double eps_inf_;
  double omega_p_;
  double gamma_d_;
  std::vector<LorentzPole> lorentz_poles_;
};

} // namespace materials
} // namespace vectorem
