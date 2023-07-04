// C/C++ header
#include <cstring>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>  // fill

// climath
#include <climath/core>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// snap
#include "thermodynamics.hpp"

static std::mutex thermo_mutex;

Real __attribute__((weak)) update_gammad(Real const& gammad, Real const q[]) {
  return gammad;
}

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
}

Thermodynamics const* Thermodynamics::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (mythermo_ == nullptr) {
    mythermo_ = Thermodynamics::fromAthenaInput(ParameterInput * pin);
  }

  return mythermo_;
}

void Thermodynamics::Destroy() {
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (Thermodynamics::mythermo_ != nullptr) {
    delete Thermodynamics::mythermo_;
    Thermodynamics::myapp_ = nullptr;
  }
}

Real Thermodynamics::GetPres(MeshBlock* pmb, int k, int j, int i) const {
  Real gm1 = pmb->peos->GetGamma() - 1;
  auto& u = pmb->phydro->u;

  Real rho = 0., fsig = 0., feps = 0.;
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += u(n, k, j, i);
    fsig += u(n, k, j, i) * cv_ratio_mass_[n];
    feps += u(n, k, j, i) / mu_ratio_[n];
  }

  Real KE = 0.5 * (sqr(u[IM1]) + sqr(u[IM2]) + sqr(u[IM3])) / rho;
  return gm1 * (u[IEN] - KE) * feps / fsig;
}

Real Thermodynamics::GetChi(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = updateGammad(pmb, k, j, i);
  auto& w = pmb->phydro->w;

  Real qsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
    feps += w(n, k, j, i) * (1. / mu_ratio_[n] - 1.);
  }

  return (gammad - 1.) / gammad * feps / qsig;
}

Real Thermodynamics::GetGamma(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = updateGammad(k, j, i);
  auto& w = pmb->phydro->w;

  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w(n, k, j, i) * (cv_ratio_mass_[n] - 1.);
    feps += w(n, k, j, i) * (1. / mu_ratio_[n] - 1.);
  }
  return 1. + (gammad - 1.) * feps / fsig;
}

Real Thermodynamics::RovRd(MeshBlock* pmb, int k, int j, int i) const {
  auto& w = pmb->phydro->w;

  Real feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) feps += w[n] * (1. / mu_ratio_[n] - 1.);
  return feps;
}

Real Thermodynamics::MoistStaticEnergy(MeshBlock* pmb, Real gz, int k, int j,
                                       int i) const {
  auto& w = pmb->phydro->w;
  auto& q = pmb->pcloud->w;

  Real temp = GetTemp(k, j, i);
  Real IE = w(IDN, k, j, i) * GetCpMass(n, k, j, i) * temp;
  Real rho_gas = w(IDN, k, j, i);
  Real rho_total = rho_gas;

  for (int n = 0; n < NCLOUD; ++n) {
    IE += rho * w(n, k, j, i) * GetCpMass(n, k, j, i) * temp;
    rho_total += rho * w(n, k, j, i);
  }

  return IE / rho_total + gz;
}

Real Thermodynamics::GetCpMass(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = updateGammad(pmb, k, j, i);
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
  return gammad / (gammad - 1.) * Rd_ * qsig;
}

Real Thermodynamics::GetCvMass(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = updateGammad(pmb, k, j, i);
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cv_ratio_mass_[n] - 1.);
  return 1. / (gammad - 1.) * Rd_ * qsig;
}

Real Thermodynamics::GetEnthalpyMass(MeshBlock* pmb, int k, int j,
                                     int i) const {
  Real gammad = updateGammad(pmb, k, j, i);
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
  return gammad / (gammad - 1.) * Rd_ * qsig * tem[0];
}

Real Thermodynamics::GetMeanMolecularWeight(MeshBlock* pmb, int k, int j,
                                            int i) const {
  auto& w = pmb->phydro->w;

  Real feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) feps += w(n, k, j, i) * (mu_ratio_[n] - 1.);
  return mu_[0] * feps;
}

Real Thermodynamics::RelativeHumidity(MeshBlock* pmb, int n, int k, int j,
                                      int i) const {
  Real dw[1 + NVAPOR];
  auto& w = pmb->phydro->w;

  Variable var;
  Real dw[1 + NVAPOR];

  pmb->pimpl->GatherPrimitive(&var, k, j, i);
  SaturationSurplus(dw, &var);
  return w(n, k, j, i) / (w(n, k, j, i) + dw[n]);
}

void Thermodynamics::SaturationSurplus(Real dw[], Variable* var) const {
  var.ConvertTo(Variable::Type::MoleFraction);

  for (int iv = 1; iv <= NVAPOR; ++iv) {
    auto rates = TryEquilibriumTP(var, iv, 0., true.);
    dw[iv] = rates[0] / var->w[iv] * var.w[iv];
  }
}

/*void Thermodynamics::UpdateTPConservingU(Real q[], Real rho, Real uhat) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real cv = 1., qtol = 1., qeps = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    uhat += beta_[n]*t3_[n]*q[n];
    qtol -= q[n];
  }
  for (int n = 1; n < NMASS; ++n) {
    cv += (cv_ratios_[n]*mu_ratios_[n] - 1.)*q[n];
    qeps += q[n]*(mu_ratios_[n] - 1.);
  }
  q[IDN] = (gamma - 1.)*uhat/cv;
  q[IPR] = rho*Rd_*q[IDN]*qtol/qeps;
}*/
