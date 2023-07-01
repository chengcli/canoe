// C/C++ header
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>  // fill

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// snap
#include "thermodynamics.hpp"

std::ostream& operator<<(std::ostream& os, Thermodynamics const& my) {
  os << "Rd [J/kg]: " << my.Rd_ << std::endl;
  os << "eps: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.mu_ratios_[i] << " ";
  os << std::endl;
  os << "cp ratios: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.cp_ratios_[i] << " ";
  os << std::endl;
  os << "cv ratios: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.cv_ratios_[i] << " ";
  os << std::endl;
  os << "beta: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.beta_[i] << " ";
  os << std::endl;
  os << "delta: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.delta_[i] << " ";
  os << std::endl;
  os << "Latent [J/kg]: ";
  for (int i = 0; i < 1 + 3 * NVAPOR; ++i) os << my.latent_[i] << " ";
  os << std::endl;
  os << "Ttriple [K]: ";
  for (int i = 0; i < NVAPOR; ++i) os << my.t3_[i] << " ";
  os << std::endl;
  os << "Ptriple [pa]: ";
  for (int i = 0; i < NVAPOR; ++i) os << my.p3_[i] << " ";
  os << std::endl;

  return os;
}

void ReadThermoProperty(Real var[], char const name[], int len, Real v0,
                        ParameterInput* pin) {
  std::stringstream msg;
  char buf[80], cstr[1024], *p;

  var[0] = v0;
  for (int n = 1; n <= NVAPOR; ++n) {
    snprintf(buf, sizeof(buf), "%s%d", name, n);
    std::string str = pin->GetString("thermodynamics", buf);
    std::snprintf(cstr, sizeof(cstr), "%s", str.c_str());
    p = std::strtok(cstr, " ,");
    int m = 0;
    while ((p != NULL) && (m++ < len)) {
      var[n + (m - 1) * NVAPOR] = std::stod(p);
      p = std::strtok(NULL, " ,");
    }
    if (m != len) {
      msg << "### FATAL ERROR in function ReadThermoProperty" << std::endl
          << "Length of '" << name << "' "
          << "doesn't equal to " << len;
      ATHENA_ERROR(msg);
    }
  }
}

Thermodynamics::Thermodynamics(MeshBlock* pmb, ParameterInput* pin)
    : pmy_block_(pmb) {
  Application::Logger app("snap");
  app->Log("Initialize Thermodynamics");

  Rd_ = pin->GetOrAddReal("thermodynamics", "Rd", 1.);

  Real gamma = pin->GetReal("hydro", "gamma");

  // Read molecular weight ratios
  ReadThermoProperty(mu_ratios_, "eps", 3, 1., pin);

  // Read cp ratios
  ReadThermoProperty(cp_ratios_, "rcp", 3, 1., pin);

  // Read beta parameter
  ReadThermoProperty(beta_, "beta", 3, 0., pin);

  // Read triple point temperature
  ReadThermoProperty(t3_, "Ttriple", 1, 0., pin);

  // Read triple point pressure
  ReadThermoProperty(p3_, "Ptriple", 1, 0., pin);

  // calculate latent heat = $\beta\frac{R_d}{\epsilon}T^r$
  for (int n = 0; n <= NVAPOR; ++n) latent_[n] = 0.;
#if NVAPOR > 0
  for (int n = 1 + NVAPOR; n < 1 + 3 * NVAPOR; ++n)
    latent_[n] = beta_[n] * Rd_ / mu_ratios_[n] * t3_[1 + (n - 1) % NVAPOR];
#endif

  // calculate delta = $(\sigma_j - \sigma_i)*\epsilon_i*\gamma/(\gamma - 1)$
  for (int n = 0; n <= NVAPOR; ++n) delta_[n] = 0.;
#if NVAPOR > 0
  for (int n = 1 + NVAPOR; n < 1 + 3 * NVAPOR; ++n)
    delta_[n] = (cp_ratios_[n] - cp_ratios_[1 + (n - 1) % NVAPOR]) *
                mu_ratios_[n] / (1. - 1. / gamma);
#endif

  // calculate cv_ratios = $\gamma\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n)
    cv_ratios_[n] = gamma * cp_ratios_[n] + (1. - gamma) / mu_ratios_[n];
  for (int n = 1 + NVAPOR; n < 1 + 3 * NVAPOR; ++n)
    cv_ratios_[n] = gamma * cp_ratios_[n];

  ftol_ = pin->GetOrAddReal("thermodynamics", "ftol", 1.0E-4);
  max_iter_ = pin->GetOrAddInteger("thermodynamics", "max_iter", 10);
}

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
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
}

void Thermodynamics::ConservedToThermodynamic(AthenaArray<Real> &q,
  AthenaArray<Real> const& u, int il, int iu, int jl, int ju, int kl, int ku)
const
{
  Real gamma = pmy_block_->peos->GetGamma();

  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i) {
        // density to molar mixing ratio
        Real rho = 0., rho_hat = 0.;
        for (int n = 0; n < NMASS; ++n) {
          rho += u(n,k,j,i);
          rho_hat += u(n,k,j,i)/mu_ratios_[n];
          q(n,k,j,i) = u(n,k,j,i)/mu_ratios_[n];
        }
        for (int n = 0; n < NMASS; ++n)
          q(n,k,j,i) /= rho_hat;

        // calculate internal energy
        Real KE = 0.5*(u(IM1,k,j,i)*u(IM1,k,j,i)
                     + u(IM2,k,j,i)*u(IM2,k,j,i)
                     + u(IM3,k,j,i)*u(IM3,k,j,i))/rho;
        Real uhat = (u(IEN,k,j,i) - KE)/(Rd_*rho_hat);

        // calculate temperature and pressure
        Real cv = 1., qtol = 1., qeps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          uhat += beta_[n]*t3_[n]*q(n,k,j,i);
          qtol -= q(n,k,j,i);
        }
        for (int n = 1; n < NMASS; ++n) {
          cv += (cv_ratios_[n]*mu_ratios_[n] - 1.)*q(n,k,j,i);
          qeps += q(n,k,j,i)*(mu_ratios_[n] - 1.);
        }
        q(IDN,k,j,i) = (gamma - 1.)*uhat/cv;
        q(IPR,k,j,i) = rho*Rd_*q(IDN,k,j,i)*qtol/qeps;
      }
}

void Thermodynamics::ThermodynamicToConserved(AthenaArray<Real> &u,
  AthenaArray<Real> const& q, int il, int iu, int jl, int ju, int kl, int ku)
const
{
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i) {
        // calculate total density
        Real qtol = 1., qeps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n)
          qtol -= q(n,k,j,i);
        for (int n = 1; n < NMASS; ++n)
          qeps += q(n,k,j,i)*(mu_ratios_[n] - 1.);
        Real rho = (q(IPR,k,j,i)*qeps)/(Rd_*q(IDN,k,j,i)*qtol);

        // molar mixing ratio to density
        Real sum = 1.;
        for (int n = 1; n < NMASS; ++n)
          sum += q(n,k,j,i)*(mu_ratios_[n] - 1.);
        for (int n = 1; n < NMASS; ++n)
          u(n,k,j,i) = rho*q(n,k,j,i)*mu_ratios_[n]/sum;
      }
}

void Thermodynamics::PolytropicIndex(AthenaArray<Real> &gm, AthenaArray<Real>
&w, int kl, int ku, int jl, int ju, int il, int iu) const
{
  Real gamma = pmy_block_->peos->GetGamma();

  // calculate local polytropic index
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real fsig = 1., feps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          fsig += w(n,k,j,i)*(cv_ratios_[n] - 1.);
          feps -= w(n,k,j,i);
        }
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += w(n,k,j,i)*(cv_ratios_[n] - 1.);
          feps += w(n,k,j,i)*(1./mu_ratios_[n] - 1.);
        }
        gm(k,j,i) = 1. + (gamma - 1.)*feps/fsig;
      }
}*/
