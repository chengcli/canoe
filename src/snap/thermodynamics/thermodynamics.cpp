// C/C++ header
#include <cstring>
#include <fstream>
#include <iomanip>
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
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// climath
#include <climath/core.h>

// snap
#include "molecules.hpp"
#include "thermodynamics.hpp"

static std::mutex thermo_mutex;

Real __attribute__((weak))
Thermodynamics::GetGammad(AirParcel const& qfrac) const {
  return gammad_ref_;
}

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
}

Thermodynamics* Thermodynamics::fromYAMLInput(YAML::Node node) {
  mythermo_ = new Thermodynamics();

  return mythermo_;
}

Thermodynamics* Thermodynamics::fromLegacyInput(ParameterInput* pin) {
  if (NCLOUD < (NPHASE - 1) * NVAPOR) {
    throw RuntimeError(
        "Thermodynamics",
        "NCLOUD < (NPHASE-1)*NVAPOR is not supported for legacy input");
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
  read_thermo_property(mythermo_->t3_.data(), "Ttriple", 1, 0., pin);

  // Read triple point pressure
  read_thermo_property(mythermo_->p3_.data(), "Ptriple", 1, 0., pin);

  mythermo_->cloud_index_set_.resize(1 + NVAPOR);
  mythermo_->svp_func1_.resize(1 + NVAPOR);

  for (int i = 1; i <= NVAPOR; ++i) {
    mythermo_->cloud_index_set_[i].resize(NPHASE - 1);
    mythermo_->svp_func1_[i].resize(NPHASE - 1);
    for (int j = 1; j < NPHASE; ++j) {
      mythermo_->cloud_index_set_[i][j - 1] = (j - 1) * NVAPOR + (i - 1);
    }
  }

  mythermo_->Rd_ = pin->GetOrAddReal("thermodynamics", "Rd", 1.);
  mythermo_->gammad_ref_ = pin->GetReal("hydro", "gamma");

  return mythermo_;
}

Thermodynamics const* Thermodynamics::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (mythermo_ == nullptr) {
    mythermo_ = new Thermodynamics();
  }

  return mythermo_;
}

Thermodynamics const* Thermodynamics::InitFromAthenaInput(ParameterInput* pin) {
  if (mythermo_ != nullptr) {
    throw RuntimeError("Thermodynamics", "Thermodynamics has been initialized");
  }

  Application::Logger app("snap");
  app->Log("Initialize Thermodynamics");

  if (pin->DoesParameterExist("thermodynamics", "thermo_file")) {
    std::string filename = pin->GetString("thermodynamics", "thermo_file");
    std::ifstream stream(filename);
    if (stream.good() == false) {
      throw RuntimeError("Thermodynamics",
                         "Cannot open thermodynamic file: " + filename);
    }

    mythermo_ = fromYAMLInput(YAML::Load(stream));
  } else {  // legacy input
    mythermo_ = fromLegacyInput(pin);
  }

  // enroll vapor functions
  mythermo_->enrollVaporFunctions();

  // alias
  auto& Rd = mythermo_->Rd_;
  auto& gammad = mythermo_->gammad_ref_;
  auto& mu_ratio = mythermo_->mu_ratio_;
  auto& inv_mu_ratio = mythermo_->inv_mu_ratio_;
  auto& mu = mythermo_->mu_;
  auto& inv_mu = mythermo_->inv_mu_;

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
  auto& cloud_index_set = mythermo_->cloud_index_set_;

  // molecular weight
  mu[0] = Constants::Rgas / Rd;

  // inverse mu ratio
  for (int n = 0; n < mu_ratio.size(); ++n) {
    inv_mu_ratio[n] = 1. / mu_ratio[n];
    mu[n] = mu[0] * mu_ratio[n];
    inv_mu[n] = 1. / mu[n];
    cp_ratio_mole[n] = cp_ratio_mass[n] * mu_ratio[n];
  }

  // calculate latent energy = $\beta\frac{R_d}{\epsilon}T^r$
  for (int n = 0; n <= NVAPOR; ++n) {
    latent_energy_mass[n] = 0.;
    latent_energy_mole[n] = 0.;
  }

  for (int i = 1; i <= NVAPOR; ++i) {
    for (int j = 0; j < NPHASE - 1; ++j) {
      int n = cloud_index_set[i][j] + 1 + NVAPOR;
      latent_energy_mass[n] = beta[n] * Rd / mu_ratio[n] * t3[i];
      latent_energy_mole[n] = latent_energy_mass[n] * mu[n];
    }
  }

  // calculate delta = $(\sigma_j - \sigma_i)*\epsilon_i*\gamma/(\gamma - 1)$
  for (int n = 0; n <= NVAPOR; ++n) delta[n] = 0.;
  for (int i = 1; i <= NVAPOR; ++i) {
    for (int j = 0; j < cloud_index_set[i].size(); ++j) {
      int n = cloud_index_set[i][j] + 1 + NVAPOR;
      delta[n] = (cp_ratio_mass[n] - cp_ratio_mass[i]) * mu_ratio[i] /
                 (1. - 1. / gammad);
    }
  }

  // set up extra clouds
  auto pindex = IndexMap::GetInstance();

  for (int n = (NPHASE - 1) * NVAPOR; n < NCLOUD; ++n) {
    Molecule mol;
    mol.LoadThermodynamicFile(pindex->GetCloudName(n) + ".thermo");

    int j = 1 + NVAPOR + n;
    mu[j] = mol.mu();
    inv_mu[j] = 1. / mu[j];

    mu_ratio[j] = mol.mu() / mu[0];
    inv_mu_ratio[j] = 1. / mu_ratio[j];

    cp_ratio_mass[j] = mol.cp() / (mol.mu() * mythermo_->GetCpMassRef(0));
    cp_ratio_mole[j] = cp_ratio_mass[j] * mu_ratio[j];

    latent_energy_mole[j] = mol.latent();
    latent_energy_mass[j] = mol.latent() / mol.mu();

    beta[j] = mol.latent() / (Constants::Rgas * Thermodynamics::RefTemp);
    delta[j] = 0.;
  }

  // finally, set up cv ratio
  // calculate cv_ratio = $\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n) {
    cv_ratio_mass[n] = gammad * cp_ratio_mass[n] + (1. - gammad) / mu_ratio[n];
    cv_ratio_mole[n] = cv_ratio_mass[n] * mu_ratio[n];
  }

  for (int n = 1 + NVAPOR; n < Size; ++n) {
    cv_ratio_mass[n] = gammad * cp_ratio_mass[n];
    cv_ratio_mole[n] = cv_ratio_mass[n] * mu_ratio[n];
  }

  // saturation adjustment parmaeters
  mythermo_->sa_max_iter_ =
      pin->GetOrAddReal("thermodynamics", "sa.max_iter", 10);
  mythermo_->sa_relax_ = pin->GetOrAddReal("thermodynamics", "sa.relax", 0.8);
  mythermo_->sa_ftol_ = pin->GetOrAddReal("thermodynamics", "sa.ftol", 1.e-4);

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
#pragma omp simd reduction(+ : rho, fsig, feps)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += u(n, k, j, i);
    fsig += u(n, k, j, i) * cv_ratio_mass_[n];
    feps += u(n, k, j, i) * inv_mu_ratio_[n];
  }

  // TODO(cli): not correct for Cubed sphere
  Real KE =
      0.5 *
      (sqr(u(IM1, k, j, i)) + sqr(u(IM2, k, j, i)) + sqr(u(IM3, k, j, i))) /
      rho;
  return gm1 * (u(IEN, k, j, i) - KE) * feps / fsig;
}

// Eq.71 in Li2019
Real Thermodynamics::GetChi(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = pmb->peos->GetGamma();
  auto& w = pmb->phydro->w;

  Real qsig = 1., feps = 1.;
#pragma omp simd reduction(+ : qsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    feps += w(n, k, j, i) * (inv_mu_ratio_[n] - 1.);
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
  }

  return (gammad - 1.) / gammad * feps / qsig;
}

// TODO(cli): check
Real Thermodynamics::GetChi(AirParcel const& qfrac) const {
  Real gammad = GetGammad(qfrac);

  Real qsig = 1., feps = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n) {
    qsig += qfrac.w[n] * (cp_ratio_mole_[n] - 1.);
  }

#pragma omp simd reduction(+ : qsig, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    feps += -qfrac.c[n];
    qsig += qfrac.c[n] * (cp_ratio_mole_[n + 1 + NVAPOR] - 1.);
  }

  return (gammad - 1.) / gammad / qsig;
}

// Eq.63 in Li2019
Real Thermodynamics::GetGamma(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = pmb->peos->GetGamma();

  AirParcel air(AirParcel::Type::MassFrac);
  pmb->pimpl->GatherFromPrimitive(&air, k, j, i);

  Real fsig = 1., feps = 1.;
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += air.w[n] * (cv_ratio_mass_[n] - 1.);
    feps += air.w[n] * (inv_mu_ratio_[n] - 1.);
  }

  for (int n = 0; n < NCLOUD; ++n) {
    fsig += air.c[n] * (cv_ratio_mass_[n + 1 + NVAPOR] - 1.);
    feps += -air.c[n];
  }

  return 1. + (gammad - 1.) * feps / fsig;
}

// Eq.94 in Li2019
Real Thermodynamics::RovRd(AirParcel const& qfrac) const {
  Real fgas = 1., feps = 1.;

#pragma omp simd reduction(+ : feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    feps += qfrac.w[n] * (mu_ratio_[n] - 1.);
  }

#pragma omp simd reduction(+ : fgas, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    fgas += -qfrac.c[n];
    feps += qfrac.c[n] * (mu_ratio_[n] - 1.);
  }

  return fgas / feps;
}

Real Thermodynamics::MoistStaticEnergy(MeshBlock* pmb, Real gz, int k, int j,
                                       int i) const {
  AirParcel var(AirParcel::Type::MassConc);
  pmb->pimpl->GatherFromPrimitive(&var, k, j, i);

  Real temp = GetTemp(pmb, k, j, i);
  Real rho = 0., LE = 0., IE = 0.;

#pragma omp simd reduction(+ : IE, rho)
  for (int n = 0; n <= NVAPOR; ++n) {
    IE += var.w[n] * GetCpMassRef(n) * temp;
    rho += var.w[n];
  }

#pragma omp simd reduction(+ : IE, LE, rho)
  for (int n = 0; n < NCLOUD; ++n) {
    IE += var.c[n] * GetCpMassRef(1 + NVAPOR + n) * temp;
    LE += -latent_energy_mass_[1 + NVAPOR + n] * var.c[n];
    rho += var.c[n];
  }

  return (IE + LE) / rho + gz;
}

// TODO(cli): check
Real Thermodynamics::GetCpMass(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = pmb->peos->GetGamma();
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
  return gammad / (gammad - 1.) * Rd_ * qsig;
}

// TODO(cli): check
Real Thermodynamics::GetCvMass(MeshBlock* pmb, int k, int j, int i) const {
  Real gammad = pmb->peos->GetGamma();
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cv_ratio_mass_[n] - 1.);
  return 1. / (gammad - 1.) * Rd_ * qsig;
}

Real Thermodynamics::GetEnthalpyMass(MeshBlock* pmb, int k, int j,
                                     int i) const {
  Real gammad = pmb->peos->GetGamma();
  auto& w = pmb->phydro->w;

  Real qsig = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n)
    qsig += w(n, k, j, i) * (cp_ratio_mass_[n] - 1.);
  return gammad / (gammad - 1.) * Rd_ * qsig * w(IDN, k, j, i);
}

Real Thermodynamics::GetMu(MeshBlock* pmb, int k, int j, int i) const {
  auto& w = pmb->phydro->w;

  Real feps = 1.;
#pragma omp simd reduction(+ : feps)
  for (int n = 1; n <= NVAPOR; ++n) feps += w(n, k, j, i) * (mu_ratio_[n] - 1.);
  return mu_[0] * feps;
}

Real Thermodynamics::RelativeHumidity(MeshBlock* pmb, int n, int k, int j,
                                      int i) const {
  AirParcel qfrac(AirParcel::Type::MoleFrac);
  pmb->pimpl->GatherFromPrimitive(&qfrac, k, j, i);
  return RelativeHumidity(qfrac, n);
}

void Thermodynamics::Extrapolate(AirParcel* qfrac, Real dzORdlnp, Method method,
                                 Real grav, Real userp) const {
  // equilibrate vapor with clouds
  for (int i = 1; i <= NVAPOR; ++i) {
    auto rates = TryEquilibriumTP_VaporCloud(*qfrac, i);

    // vapor condensation rate
    qfrac->w[i] += rates[0];

    // cloud concentration rates
    for (int j = 1; j < rates.size(); ++j)
      qfrac->c[cloud_index_set_[i][j - 1]] += rates[j];
  }

  // RK4 integration
#ifdef HYDROSTATIC
  rk4IntegrateLnp(qfrac, dzORdlnp, method, userp);
#else
  rk4IntegrateZ(qfrac, dzORdlnp, method, grav, userp);
#endif
}

Real Thermodynamics::GetLatentHeatMole(int i, std::vector<Real> const& rates,
                                       Real temp) const {
  if (std::abs(rates[0]) < 1.E-8) return 0.;

  Real heat = 0.;
  for (int j = 1; j < rates.size(); ++j) {
    int n = cloud_index_set_[i][j - 1] + 1 + NVAPOR;
    heat += rates[j] * GetLatentEnergyMole(n, temp);
  }

  return heat / std::abs(rates[0]);
}

// Eq.4.5.11 in Emanuel (1994)
Real Thermodynamics::EquivalentPotentialTemp(MeshBlock* pmb, Real p0, int v,
                                             int k, int j, int i) const {
#if (NVAPOR > 0)
  AirParcel air_mass(AirParcel::Type::MassFrac);
  pmb->pimpl->GatherFromPrimitive(&air_mass, k, j, i);

  AirParcel air_mole(AirParcel::Type::MoleFrac);
  pmb->pimpl->GatherFromPrimitive(&air_mole, k, j, i);

  // get dry air mixing ratio
  Real xg = 1., qd = 1.;
#pragma omp simd reduction(+ : xg, qd)
  for (int n = 0; n < NCLOUD; ++n) {
    xg += -air_mole.c[n];
    qd += -air_mass.c[n];
  }

  Real xd = xg;
#pragma omp simd reduction(+ : xd, qd)
  for (int n = 1; n <= NVAPOR; ++n) {
    xd += -air_mole.w[n];
    qd += -air_mass.w[n];
  }

  Real temp = air_mole.w[IDN];
  Real pres = air_mole.w[IPR];

  Real qv = air_mass.w[v];
  Real qc = air_mass.c[cloud_index_set_[v][0]];
  Real qt = qv + qc;

  Real Rd = Rd_;
  Real Rv = Rd_ / mu_ratio_[v];

  Real cpd = GetCpMassRef(0);
  Real cl = GetCpMassRef(cloud_index_set_[v][0] + 1 + NVAPOR);

  std::vector<Real> rates{-0.01, 0.01};
  Real lv = GetLatentHeatMass(v, rates, temp);

  Real rh = RelativeHumidity(pmb, v, k, j, i);
  Real pd = xd / xg * pres;
  Real cpt = cpd * qd + cl * qt;

  return temp * pow(p0 / pd, Rd * qd / cpt) * pow(rh, -Rv * qv / cpt) *
         exp(lv * qv / (cpt * temp));
#else
  return PotentialTemp(pmb, p0, k, j, i);
#endif
}

void Thermodynamics::updateTPConservingU(AirParcel* qfrac, Real rmole,
                                         Real umole) const {
  Real cvd = Constants::Rgas / (GetGammad(*qfrac) - 1.);
  Real fsig = 1., qgas = 1.;

  for (int i = 1; i <= NVAPOR; ++i) {
    // vapor
    fsig += (cv_ratio_mole_[i] - 1.) * qfrac->w[i];

    // clouds
    for (auto j : cloud_index_set_[i]) {
      int n = j + 1 + NVAPOR;
      Real qc = qfrac->c[j];

      fsig += (cv_ratio_mole_[n] - 1.) * qc;
      umole += latent_energy_mole_[n] * qc;
    }
  }

  // clouds
#pragma omp simd reduction(+ : qgas)
  for (int n = 0; n < NCLOUD; ++n) qgas += -qfrac->c[n];

  qfrac->w[IDN] = umole / (cvd * fsig);
  qfrac->w[IPR] = rmole * qgas * Constants::Rgas * qfrac->w[IDN];
}

Real Thermodynamics::getInternalEnergyMole(AirParcel const& qfrac) const {
  Real cvd = Constants::Rgas / (GetGammad(qfrac) - 1.);
  Real fsig = 1., LE = 0.;

  for (int i = 1; i <= NVAPOR; ++i) {
    // vapor
    fsig += (cv_ratio_mole_[i] - 1.) * qfrac.w[i];

    // clouds
    for (auto j : cloud_index_set_[i]) {
      int n = j + 1 + NVAPOR;
      Real qc = qfrac.c[j];

      fsig += (cv_ratio_mole_[n] - 1.) * qc;
      LE += latent_energy_mole_[n] * qc;
    }
  }

  return cvd * qfrac.w[IDN] * fsig - LE;
}

Real Thermodynamics::getDensityMole(AirParcel const& qfrac) const {
  Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
  for (int n = 0; n < NCLOUD; ++n) qgas += -qfrac.c[n];
  return qfrac.w[IPR] / (Constants::Rgas * qfrac.w[IDN] * qgas);
}

void Thermodynamics::setTotalEquivalentVapor(AirParcel* qfrac) const {
  // vpaor <=> cloud
  for (int i = 1; i <= NVAPOR; ++i)
    for (auto& j : cloud_index_set_[i]) {
      qfrac->w[i] += qfrac->c[j];
      qfrac->c[j] = 0.;
    }

  // vapor + vapor <=> cloud
  for (auto const& [ij, info] : cloud_reaction_map_) {
    auto& indx = info.first;
    auto& stoi = info.second;

    qfrac->w[indx[0]] += stoi[0] / stoi[2] * qfrac->c[indx[2]];
    qfrac->c[indx[2]] = 0.;
    qfrac->w[indx[1]] += stoi[1] / stoi[2] * qfrac->c[indx[2]];
    qfrac->c[indx[2]] = 0.;
  }
}

Thermodynamics* Thermodynamics::mythermo_ = nullptr;
