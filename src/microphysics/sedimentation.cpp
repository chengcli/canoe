// athena
#include <athena/athena.hpp>

// shared
#include <constants.hpp>

// microphysics
#include "sedimentation.hpp"

void SedimentationImpl::reset() {
  if (options.radius().size() != options.density().size()) {
    throw std::runtime_error(
        "Sedimentation: radius and density must have the same size");
  }

  radius = register_parameter("radius", torch::tensor(options.radius()));
  density = register_parameter("density", torch::tensor(options.density()));
}

torch::Tensor SedimentationImpl::forward(torch::Tensor hydro_w) {
  const auto d = options.a_diameter();
  const auto epsilon_LJ = options.a_epsilon_LJ();
  const auto m = options.a_mass();

  auto tem = shared_->at("temperature").get();

  auto eta = (5.0 / 16.0) * torch::sqrt(M_PI * m * Constants::kBoltz * tem) *
             torch::pow(Constants::kBoltz * tem / epsilon_LJ, 0.16) /
             (M_PI * d * d * 1.22);

  // Calculate mean free path, lambda
  auto lambda =
      (eta * std::sqrt(M_PI * Constants::kBoltz * Constants::kBoltz)) /
      (hydro_w[IPR] * std::sqrt(2.0 * m));

  // Calculate Knudsen number, Kn
  auto Kn = lambda / radius;

  // Calculate Cunningham slip factor, beta
  auto beta = 1.0 + Kn * (1.256 + 0.4 * torch::exp(-1.1 / Kn));

  // Calculate vsed
  auto vel =
      beta *
      (2.0 * radius * radius * options.gravity() * (density - hydro_w[IDN])) /
      (9.0 * eta);
  //          std::cout << "vel: " << vel << " pressure:" << P
  //                    << "rho_gas: " << rho_gas << " T:" << T <<
  //                    "beta: " << beta
  //                    << " eta:" << eta << " rho_gas: " << rho_gas
  //                    << " lambda:" << lambda << std::endl;

  return torch::clamp(vel, c10::nullopt, options.upper_limit());
}
