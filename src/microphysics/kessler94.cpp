// C/C++
#include <string>

// external
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>
#include <Eigen/Dense>

// Athena++
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/parameter_map.hpp>

// canoe
#include <configure.hpp>
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// microphysics
#include "chemistry_solver.hpp"
#include "microphysical_scheme.hpp"

Kessler94::Kessler94(std::string name, YAML::Node &node)
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

Kessler94::~Kessler94() {
  Application::Logger app("microphysics");
  app->Log("Destroy Kessler94 Scheme");
}

void Kessler94::AssembleReactionMatrix(Real *rate, Real **jac,
                                       AirParcel const &air, Real time) {
  int iv = species_index_[0];
  int ic = species_index_[1];
  int ip = species_index_[2];

  for (int i : {iv, ic, ip}) {
    rate[i] = 0.;
    for (int j : {iv, ic, ip}) {
      jac[i][j] = 0.;
    }
  }

  Real k1 = params_["autoconversion"];
  Real k2 = params_["accretion"];
  Real k3 = params_["evaporation"];

  auto pthermo = Thermodynamics::GetInstance();

  /*Real svp = pmy_block_->pthermo->SatVaporPressure(air.w[IDN], ic);
  Real qtol = 1.;  // total gas mols
  for (int n = 0; n < NCLOUD; ++n) qtol -= air.c[n];
  Real dqH2O = air.w[iv] - svp / air.w[IPR] * qtol;  // q - qs*/

  Real rh = pthermo->RelativeHumidity(air, iv);
  Real factor = 1. - rh;

  if (rh < 1.) {  // evaporation
    rate[iv] += -k3 * air.w[ip] * factor;
    rate[ip] += k3 * air.w[ip] * factor;
    jac[iv][iv] += -k3 * air.w[ip];
    jac[iv][ip] += -k3 * factor;
    jac[ip][iv] += k3 * air.w[ip];
    jac[ip][ip] += k3 * factor;
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

void Kessler94::EvolveOneStep(AirParcel *air, Real time, Real dt) {
  MapVector rate(&rate_[0]);
  MapMatrix jacobian(&jacobian_[0][0]);

  auto sol = solver_.solveBDF1<RealVector>(rate, jacobian, dt);
  for (int n = 0; n < Base::Size; ++n) air->w[species_index_[n]] += sol(n);
}
