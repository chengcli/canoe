// C/C++ header
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>  // fill

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>
#include <impl.hpp>
#include <variable.hpp>

// climath
#include <climath/core.h>

// snap
#include "thermodynamics.hpp"

static std::mutex thermo_mutex;

Real __attribute__((weak)) update_gammad(Real gammad, Real const q[]) {
  return gammad;
}

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
}

Thermodynamics::Thermodynamics(YAML::Node& node) {}

Thermodynamics const* Thermodynamics::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (mythermo_ == nullptr) {
    mythermo_ = new Thermodynamics();
  }

  return mythermo_;
}

Thermodynamics const* Thermodynamics::InitFromAthenaInput(ParameterInput* pin) {
  Application::Logger app("snap");
  app->Log("Initialize Thermodynamics");

  if (pin->DoesParameterExist("thermodynamics", "control_file") == false) {
    std::string filename = pin->GetString("thermodynamics", "control_file");
    std::ifstream stream(filename);
    if (stream.good() == false) {
      throw RuntimeError("Thermodynamics",
                         "Cannot open thermodynamic file: " + filename);
    }
    YAML::Node node = YAML::Load(stream);
    mythermo_ = new Thermodynamics(node);
  } else {  // legacy input
    if (NCLOUD != 2 * NVAPOR) {
      throw RuntimeError(
          "Thermodynamics",
          "NCLOUD != 2*NVAPOR is not supported for legacy input");
    }

    mythermo_ = new Thermodynamics();

    // Read molecular weight ratios
    read_thermo_property(mythermo_->mu_ratio_.data(), "eps", NPHASE, 1., pin);

    // Read cp ratios
    read_thermo_property(mythermo_->cp_ratio_mass_.data(), "rcp", NPHASE, 1.,
                         pin);

    // Read beta parameter
    read_thermo_property(mythermo_->beta_.data(), "beta", NPHASE, 0., pin);

    // Read triple point temperature
    read_thermo_property(mythermo_->t3_.data(), "Ttriple", NPHASE, 0., pin);

    // Read triple point pressure
    read_thermo_property(mythermo_->p3_.data(), "Ptriple", NPHASE, 0., pin);

    mythermo_->cloud_index_set_.resize(NVAPOR);

    for (int i = 1; i < NVAPOR; ++i) {
      mythermo_->cloud_index_set_[i].resize(NPHASE - 1);
      for (int j = 1; j < NPHASE; ++j) {
        mythermo_->cloud_index_set_[i][j - 1] = 1 + j * NVAPOR + i - 1;
      }
    }
  }

  mythermo_->Rd_ = pin->GetOrAddReal("thermodynamics", "Rd", 1.);
  mythermo_->gammad_ = pin->GetReal("hydro", "gamma");

  // alias
  auto& Rd = mythermo_->Rd_;
  auto& gammad = mythermo_->Rd_;
  auto& mu_ratio = mythermo_->mu_ratio_;

  auto& cp_ratio_mole = mythermo_->cp_ratio_mole_;
  auto& cp_ratio_mass = mythermo_->cp_ratio_mass_;
  auto& cv_ratio_mole = mythermo_->cv_ratio_mole_;
  auto& cv_ratio_mass = mythermo_->cv_ratio_mass_;

  auto& latent_energy_mass = mythermo_->latent_energy_mass_;
  auto& latent_energy_mole = mythermo_->latent_energy_mole_;

  auto& beta = mythermo_->beta_;
  auto& delta = mythermo_->delta_;
  auto& t3 = mythermo_->t3_;
  auto& p3 = mythermo_->p3_;

  // calculate latent energy = $\beta\frac{R_d}{\epsilon}T^r$
  for (int n = 0; n <= NVAPOR; ++n) latent_energy_mass[n] = 0.;
  for (int n = 1 + NVAPOR; n < Size; ++n)
    latent_energy_mass[n] = beta[n] * Rd / mu_ratio[n] * t3[n];

  // calculate delta = $(\sigma_j - \sigma_i)*\epsilon_i*\gamma/(\gamma - 1)$
  for (int n = 0; n <= NVAPOR; ++n) delta[n] = 0.;
  for (int n = 1 + NVAPOR; n < Size; ++n)
    delta[n] = (cp_ratio_mass[n] - cp_ratio_mass[1 + (n - 1) % NVAPOR]) *
               mu_ratio[n] / (1. - 1. / gammad);

  // calculate cv_ratio = $\gamma\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n)
    cv_ratio_mass[n] = gammad * cp_ratio_mass[n] + (1. - gammad) / mu_ratio[n];
  for (int n = 1 + NVAPOR; n < 1 + 3 * NVAPOR; ++n)
    cv_ratio_mass[n] = gammad * cp_ratio_mass[n];

  return mythermo_;
}

void Thermodynamics::Destroy() {
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (Thermodynamics::mythermo_ != nullptr) {
    delete Thermodynamics::mythermo_;
    Thermodynamics::mythermo_ = nullptr;
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

  Real KE =
      0.5 *
      (sqr(u(IM1, k, j, i)) + sqr(u(IM2, k, j, i)) + sqr(u(IM3, k, j, i))) /
      rho;
  return gm1 * (u(IEN, k, j, i) - KE) * feps / fsig;
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

// TODO(cli): check
Real Thermodynamics::GetChi(Variable const& qfrac) const {
  Real gammad = updateGammad(qfrac);

  Real qsig = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    qsig += qfrac.w[n] * (cp_ratio_mole_[n] - 1.);
  }

  return (gammad - 1.) / gammad / qsig;
}

Real Thermodynamics::GetGamma(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = updateGammad(pmb, k, j, i);
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
  for (int n = 1; n <= NVAPOR; ++n)
    feps += w(n, k, j, i) * (1. / mu_ratio_[n] - 1.);
  return feps;
}

Real Thermodynamics::MoistStaticEnergy(MeshBlock* pmb, Real gz, int k, int j,
                                       int i) const {
  auto& w = pmb->phydro->w;
  auto& c = pmb->pimpl->pcloud->w;

  Real temp = GetTemp(pmb, k, j, i);
  Real IE = w(IDN, k, j, i) * GetCpMass(pmb, k, j, i) * temp;
  Real rho_gas = w(IDN, k, j, i);
  Real rho_total = rho_gas;

  for (int n = 0; n < NCLOUD; ++n) {
    IE += rho_gas * c(n, k, j, i) * GetCpMass(n) * temp;
    rho_total += rho_gas * c(n, k, j, i);
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
  return gammad / (gammad - 1.) * Rd_ * qsig * w(IDN, k, j, i);
}

Real Thermodynamics::GetMu(MeshBlock* pmb, int k, int j, int i) const {
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

  pmb->pimpl->GatherPrimitive(&var, k, j, i);
  getSaturationSurplus(dw, var);
  return w(n, k, j, i) / (w(n, k, j, i) + dw[n]);
}

void Thermodynamics::Extrapolate(Variable* qfrac, Real dzORdlnp, Method method,
                                 Real grav, Real userp) const {
  Real latent[1 + NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n) latent[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    auto rates = TryEquilibriumTP(*qfrac, iv);

    // saturation indicator
    latent[iv] = getLatentHeat(rates, iv);

    // vapor condensation rate
    qfrac->w[iv] += rates[0];

    // cloud concentration rates
    for (int n = 1; n < rates.size(); ++n)
      qfrac->c[cloud_index_set_[iv][n - 1]] += rates[n];
  }

  // RK4 integration
#ifdef HYDROSTATIC
  rk4IntegrateLnp(qfrac, latent, dzORdlnp, method, userp);
#else
  rk4IntegrateZ(qfrac, latent, dzORdlnp, method, grav, userp);
#endif
}

void Thermodynamics::getSaturationSurplus(Real dw[], Variable& var) const {
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    dw[iv] = var.w[iv];
  }

  var.ConvertToMoleFraction();
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    auto rates = TryEquilibriumTP(var, iv, 0., true);
    dw[iv] *= rates[0] / var.w[iv];
  }
}

Thermodynamics* Thermodynamics::mythermo_ = nullptr;

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
