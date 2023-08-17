// Athena++
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>

// snap
#include <index_map.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

#include "chemistry_solver.hpp"
#include "microphysical_scheme.hpp"

void Kessler94::Kessler94(YAML::Node &node) : MicrophysicalScheme<3>(node) {
  Application::Logger app("microphysics");
  app->Log("Initialize Kessler94 Scheme");

  Debugger *pdbg = pchem->pmy_block->pdebug;
  std::stringstream &msg = pdbg->msg;

  myname = name;
  particle_name = pin->GetString("chemistry", name + ".link_particle");
  msg << "- particle " << particle_name << " linked to " << name << " chemistry"
      << std::endl;

  coeffs_["condensation"] = pin->GetReal("chemistry", name + ".condensation");
  coeffs_["autoconversion"] =
      pin->GetReal("chemistry", name + ".autoconversion");
  coeffs_["accretion"] = pin->GetReal("chemistry", name + ".accretion");
  coeffs_["evaporation"] = pin->GetReal("chemistry", name + ".evaporation");

  species_index_[0] = pin->GetInteger("chemistry", name + ".link_vapor");
  species_index_[1] = NHYDRO;
  species_index_[2] = NHYDRO + 1;
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

  Real k1 = coeffs_["autoconversion"];
  Real k2 = coeffs_["accretion"];
  Real k3 = coeffs_["evaporation"];

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

void Kessler94::EvolveOneStep(AirParcel *air, Real time, Real dt) {
  Eigen::Map<Eigen::Matrix<Real, Size, 1>> rate(&rate_[0]);
  Eigen::Map<Eigen::Matrix<Real, Size, Size>> jacobian(&jacobian_[0][0]);

  auto sol = solve_.solveBDF1(rate, jac, dt);

  for (int n = 0; n < Size; ++n) air->w[species_index_[n]] += sol(n);
}
