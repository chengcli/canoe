// external
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// microphysics
#include "microphysical_schemes.hpp"

Kessler94::Kessler94(std::string name, YAML::Node const &node)
    : MicrophysicalScheme<3>(name, node) {
  Application::Logger app("microphysics");
  app->Log("Initialize Kessler94 Scheme for " + name);
}

Kessler94::~Kessler94() {
  Application::Logger app("microphysics");
  app->Log("Destroy Kessler94 Scheme for " + GetName());
}

void Kessler94::AssembleReactionMatrix(AirParcel const &air, Real time) {
  auto pthermo = Thermodynamics::GetInstance();

  // get indices
  int iv = mySpeciesId(0);
  int ic = mySpeciesId(1);
  int ip = mySpeciesId(2);

  // get parameters
  Real k1 = GetPar<Real>("autoconversion");
  Real k2 = GetPar<Real>("accretion");
  Real k3 = GetPar<Real>("evaporation");

  // calculate saturation deficit (negative means sub-saturation)
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
  auto pthermo = Thermodynamics::GetInstance();

  // auto sol = solver_.solveBDF1<Base::RealVector>(rate_, jacb_, dt);
  // auto sol = solver_.solveTRBDF2<Base::RealVector>(rate_, jacb_, dt);
  auto sol = solver_.solveTRBDF2Blend<Base::RealVector>(
      rate_, jacb_, dt, air->w, GetMySpeciesIndices().data());

  //! \todo check this
  // 0 is a special buffer place for cloud in equilibrium with vapor at the same
  // temperature
  int jbuf = pthermo->GetCloudIndex(mySpeciesId(0), 0);

  air->c[jbuf] += sol(0);
  for (int n = 1; n < Size; ++n) air->w[mySpeciesId(n)] += sol(n);
}

void Kessler94::SetVsedFromConserved(AthenaArray<Real> vsed[3],
                                     Hydro const *phydro, int kl, int ku,
                                     int jl, int ju, int il, int iu) {
  int ip = myCloudId(2);
  Real vel = GetPar<Real>("sedimentation");

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        vsed[0](ip, k, j, i) = vel;
      }
}
