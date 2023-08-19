// C/C++
#include <string>

// external
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

// Athena++
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/parameter_map.hpp>

// snap
#include <index_map.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// microphysics
#include "chemistry_solver.hpp"
#include "microphysical_scheme.hpp"

void Kessler94::Kessler94(std::string name, YAML::Node &node)
    : MicrophysicalScheme<3>(name, node) {
  Application::Logger app("microphysics");
  app->Log("Initialize Kessler94 Scheme");

  if (node["parameters"]) params_ = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  auto pindex = IndexMap::GetInstance();

  for (int i = 0; i < Base::Size; ++i) {
    species_index_[i] = pindex->GetSpeciesId(species[i]);
  }
}

~Kessler94::Kessler94() {
  Application::Logger app("microphysics");
  app->Log("Destroy Kessler94 Scheme");
}

void Kessler94::AssembleReactionMatrix(Real *rate, Real **jac,
                                       AirParcel const &air, Real time) {
  int iv = species_index_[0];
  int ic = species_index_[1];
  int ip = species_index_[2];

  Real k1 = param_["autoconversion"];
  Real k2 = param_["accretion"];
  Real k3 = param_["evaporation"];

  Real svp = pmy_block_->pthermo->SatVaporPressure(air.w[IDN], ic);
  Real qtol = 1.;  // total gas mols
  for (int n = 1 + NVAPOR; n < ITR; ++n) qtol -= air.w[n];
  Real dqH2O = air.w[iv] - svp / air.w[IPR] * qtol;  // q - qs

  if (dqH2O < 0.) {  // evaporation
    rate[iv] += -k3 * air.w[ip] * dqH2O;
    rate[ip] += k3 * air.w[ip] * dqH2O;
    jac[iv][iv] += -k3 * air.w[ip];
    jac[iv][ip] += -k3 * dqH2O;
    jac[ip][iv] += k3 * air.w[ip];
    jac[ip][ip] += k3 * dqH2O;
  }

  // autoconversion
  rate[ic] += -k1 * air.w[ic];
  rate[ip] += k1 * air.w[ic];
  jac[ic][ic] += -k1;
  jac[ip][ic] += k1;

  // accretion
  rate[ic] += -k2 * air.w[ic] * air.w[ip];
  rate[ip] += k2 * air.w[ic] * air.w[ip];
  jac[ic][ic] += -k2 * air.w[ip];
  jac[ic][ip] += -k2 * air.w[ic];
  jac[ip][ic] += k2 * air.w[ip];
  jac[ip][ip] += k2 * air.w[ic];
}

using Vector = Eigen::Map<Eigen::Matrix<Real, Base::Size, 1>>;
using Matrix = Eigen::Map<Eigen::Matrix<Real, Base::Size, Base::Size>>;

void Kessler94::EvolveOneStep(AirParcel *air, Real time, Real dt) {
  Vector rate(&rate_[0]);
  Matrix jacobian(&jacobian_[0][0]);

  auto sol = solve_.solveBDF1(rate, jac, dt);

  for (int n = 0; n < Base::Size; ++n) air->w[species_index_[n]] += sol(n);
}
