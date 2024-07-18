// C/C++ header
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>  // fill

// cantera
#include <cantera/kinetics.h>
#include <cantera/thermo.h>

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
#include "atm_thermodynamics.hpp"
#include "molecules.hpp"

static std::mutex thermo_mutex;

const std::string Thermodynamics::input_key = "thermodynamics_config";

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
}

Thermodynamics* Thermodynamics::fromYAMLInput(std::string const& fname) {
  mythermo_ = new Thermodynamics();

  auto thermo = Cantera::newThermo(fname, "gas");

  if (thermo->nSpecies() != 1 + NVAPOR + NCLOUD) {
    throw RuntimeError("Thermodynamics",
                       "Number of species does not match the input file");
  }

  auto& kinetics = mythermo_->kinetics_;
  kinetics = Cantera::newKinetics({thermo}, fname);

  // --------- vapor + cloud thermo ---------
  std::vector<Real> mu(thermo->nSpecies());
  std::vector<Real> cp_mole(thermo->nSpecies());

  // g/mol
  thermo->getMolecularWeights(mu.data());

  // J/kmol/K
  thermo->getPartialMolarCp(cp_mole.data());

  mythermo_->Rd_ = Cantera::GasConstant / mu[0];
  mythermo_->gammad_ref_ = cp_mole[0] / (cp_mole[0] - Cantera::GasConstant);

  // ---------- dimensionless properties ----------

  // cp ratios
  for (size_t i = 0; i < cp_mole.size(); ++i) {
    mythermo_->cp_ratio_mole_[i] = cp_mole[i] / cp_mole[0];
  }

  // molecular weight ratios
  for (size_t i = 0; i < mu.size(); ++i) {
    mythermo_->mu_ratio_[i] = mu[i] / mu[0];
  }

  // beta parameter (dimensionless internal energy)
  for (size_t i = 0; i < thermo->nSpecies(); ++i) {
    auto& tp = thermo->species(i)->thermo;
    size_t n;
    int type;
    Real tlow, thigh, pref;
    std::vector<Real> coeffs(tp->nCoeffs());
    tp->reportParameters(n, type, tlow, thigh, pref, coeffs.data());

    Real t0 = coeffs[0];
    Real h0 = coeffs[1];
    Real s0 = coeffs[2];
    Real cp0 = coeffs[3];
    mythermo_->beta_[i] = h0 / (Cantera::GasConstant * t0);
  }

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

  mythermo_->svp_func1_.resize(1 + NVAPOR);
  mythermo_->cloud_index_set_.resize(1 + NVAPOR);

  for (int i = 0; i <= NVAPOR; ++i) {
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

Thermodynamics const* Thermodynamics::InitFromYAMLInput(
    std::string const& fname) {
  if (mythermo_ != nullptr) {
    throw RuntimeError("Thermodynamics", "Thermodynamics has been initialized");
  }

  Application::Logger app("snap");
  app->Log("Initialize Thermodynamics");

  return fromYAMLInput(fname);
}

Thermodynamics const* Thermodynamics::InitFromAthenaInput(ParameterInput* pin) {
  if (mythermo_ != nullptr) {
    throw RuntimeError("Thermodynamics", "Thermodynamics has been initialized");
  }

  Application::Logger app("snap");
  app->Log("Initialize Thermodynamics");

  if (pin->DoesParameterExist("thermodynamics", input_key)) {
    std::string filename = pin->GetString("thermodynamics", input_key);
    mythermo_ = fromYAMLInput(filename);
  } else {  // legacy input
    mythermo_ = fromLegacyInput(pin);
  }

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
    cp_ratio_mass[n] = cp_ratio_mole[n] * inv_mu_ratio[n];
  }

  /* calculate latent energy = $\beta\frac{R_d}{\epsilon}T^r$
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
  }*/

  // finally, set up cv ratio
  // calculate cv_ratio = $\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n) {
    cv_ratio_mass[n] =
        gammad * cp_ratio_mass[n] + (1. - gammad) * inv_mu_ratio[n];
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

// Eq.71 in Li2019
Real Thermodynamics::GetChi(MeshBlock const* pmb, int k, int j, int i) const {
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

Real Thermodynamics::MoistStaticEnergy(MeshBlock const* pmb, Real gz, int k,
                                       int j, int i) const {
  auto&& air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);
  air.ToMassConcentration();

  Real temp = GetTemp(pmb, k, j, i);
  Real rho = 0., LE = 0., IE = 0.;

#pragma omp simd reduction(+ : IE, rho)
  for (int n = 0; n <= NVAPOR; ++n) {
    IE += air.w[n] * GetCpMassRef(n) * temp;
    rho += air.w[n];
  }

#pragma omp simd reduction(+ : IE, LE, rho)
  for (int n = 0; n < NCLOUD; ++n) {
    IE += air.c[n] * GetCpMassRef(1 + NVAPOR + n) * temp;
    LE += -latent_energy_mass_[1 + NVAPOR + n] * air.c[n];
    rho += air.c[n];
  }

  return (IE + LE) / rho + gz;
}

Real Thermodynamics::RelativeHumidity(MeshBlock const* pmb, int n, int k, int j,
                                      int i) const {
  auto&& air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);
  air.ToMoleFraction();
  return get_relative_humidity(air, n);
}

template <typename T>
void Thermodynamics::Extrapolate(T w, Real dzORdlnp, std::string method,
                                 Real grav, Real userp) const {
  auto& thermo = kinetics_->thermo();

  double pres = w[IPR];
  thermo.setMassFractionsPartial(w);
  thermo.setDensity(w[IDN]);
  thermo.setPressure(w[IPR]);

  EquilibrateTP();

  Real qfrac[NHYDRO];
  thermo.getMoleFractionsPartial(qfrac);
  qfrac[IDN] = thermo.temperature();
  qfrac[IPR] = thermo.pressure();

  // RK4 integration
#ifdef HYDROSTATIC
  rk4_integrate_lnp(qfrac, dzORdlnp, method, userp);
#else
  rk4_integrate_z(qfrac, dzORdlnp, method, grav, userp);
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

Thermodynamics* Thermodynamics::mythermo_ = nullptr;
