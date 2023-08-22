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

Kessler94::Kessler94(std::string name, YAML::Node const &node)
    : MicrophysicalScheme<3>(name, node) {
  Application::Logger app("microphysics");
  app->Log("Initialize Kessler94 Scheme for " + name);

  if (node["parameters"]) params_ = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  auto pindex = IndexMap::GetInstance();

  for (int i = 0; i < Size; ++i) {
    species_index_[i] = pindex->GetSpeciesId(species[i]);
  }
}

Kessler94::~Kessler94() {
  Application::Logger app("microphysics");
  app->Log("Destroy Kessler94 Scheme for " + name_);
}

void Kessler94::AssembleReactionMatrix(AirParcel const &air, Real time) {
  auto pthermo = Thermodynamics::GetInstance();

  // get indices
  int iv = species_index_[0];
  int ic = species_index_[1];
  int ip = species_index_[2];

  // get parameters
  Real k1 = params_["autoconversion"];
  Real k2 = params_["accretion"];
  Real k3 = params_["evaporation"];

  // calculate saturatin deficit (negative means sub-saturation)
  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iv, 0., true);
  Real dqv = -rates[0];

  // assemble matrix
  rate_.setZero();
  jacb_.setZero();

  if (dqv < 0.) {  // evaporation
    rate_(0) += -k3 * air.w[ip] * dqv;
    rate_(2) += k3 * air.w[ip] * dqv;
    jacb_(0, 0) += -k3 * air.w[ip];
    jacb_(0, 2) += -k3 * dqv;
    jacb_(2, 0) += k3 * air.w[ip];
    jacb_(2, 2) += k3 * dqv;
  }

  // autoconversion
  rate_(1) += -k1 * air.w[ic];
  rate_(2) += k1 * air.w[ic];
  jacb_(1, 1) += -k1;
  jacb_(2, 1) += k1;

  // accretion
  rate_(1) += -k2 * air.w[ic] * air.w[ip];
  rate_(2) += k2 * air.w[ic] * air.w[ip];
  jacb_(1, 1) += -k2 * air.w[ip];
  jacb_(1, 2) += -k2 * air.w[ic];
  jacb_(2, 1) += k2 * air.w[ip];
  jacb_(2, 2) += k2 * air.w[ic];
}

void Kessler94::EvolveOneStep(AirParcel *air, Real time, Real dt) {
  // auto sol = solver_.solveBDF1<Base::RealVector>(rate_, jacb_, dt);
  // auto sol = solver_.solveTRBDF2<Base::RealVector>(rate_, jacb_, dt);
  auto sol = solver_.solveTRBDF2Blend<Base::RealVector>(
      rate_, jacb_, dt, air->w, species_index_.data());

  for (int n = 0; n < Size; ++n) air->w[species_index_[n]] += sol(n);
}

void Kessler94::SetSedimentationVelocity(AthenaArray<Real> &vsed, int k, int j,
                                         int il, int iu) {
  int ip = species_index_[2] - NHYDRO;
  Real vel = params_["sedimentation"];

  for (int i = il; i <= iu; ++i) vsed(ip, k, j, il, iu) = vel;
}
