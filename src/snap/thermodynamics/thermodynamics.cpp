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
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>
#include <impl.hpp>

// snap
#include "thermodynamics.hpp"

static std::mutex thermo_mutex;

const std::string Thermodynamics::input_key = "thermodynamics_config";

Thermodynamics::~Thermodynamics() {
  Application::Logger app("snap");
  app->Log("Destroy Thermodynamics");
}

void Thermodynamics::UpdateThermoProperties() {
  auto& thermo = kinetics_->thermo();

  // --------- vapor + cloud thermo ---------
  std::vector<Real> mu(thermo.nSpecies());
  std::vector<Real> cp_mole(thermo.nSpecies());

  // g/mol
  thermo.getMolecularWeights(mu.data());

  // J/kmol/K
  thermo.getPartialMolarCp(cp_mole.data());

  Rd_ = Cantera::GasConstant / mu[0];
  gammad_ = cp_mole[0] / (cp_mole[0] - Cantera::GasConstant);

  // ---------- dimensionless properties ----------
  // molecular weight ratios
  for (size_t i = 0; i < mu.size(); ++i) {
    inv_mu_ratio_[i] = mu[0] / mu[i];
  }

  // cp ratios
  for (size_t i = 0; i < cp_mole.size(); ++i) {
    cp_ratio_[i] = cp_mole[i] / cp_mole[0] * inv_mu_ratio_[i];
  }

  // calculate cv_ratio = $\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n) {
    cv_ratio_[n] = gammad_ * cp_ratio_[n] + (1. - gammad_) * inv_mu_ratio_[n];
  }

  for (int n = 1 + NVAPOR; n < cv_ratio_.size(); ++n) {
    cv_ratio_[n] = gammad_ * cp_ratio_[n];
  }
}

Thermodynamics* Thermodynamics::fromYAMLInput(std::string const& fname) {
  mythermo_ = new Thermodynamics();

  auto atm = Cantera::newThermo(fname, "atm");
  if (atm->nSpecies() != 1 + NVAPOR + NCLOUD + NPRECIP) {
    throw RuntimeError("Thermodynamics",
                       "Number of species does not match the input file");
  }
  // atm->setNDim(2);

  /*auto precip = Cantera::newThermo(fname, "precip");
  if (precip->nSpecies() != NPRECIP) {
    throw RuntimeError("Thermodynamics",
                       "Number of species does not match the input file");
  }

  std::cout << "atm dim = " << atm->nDim() << std::endl;
  std::cout << "precip dim = " << precip->nDim() << std::endl;*/

  auto& kinetics = mythermo_->kinetics_;
  kinetics = std::static_pointer_cast<Cantera::Condensation>(
      Cantera::newKinetics({atm}, fname));

  // finalize setup the thermo manager for clouds
  atm->updateFromKinetics(*kinetics);
  // precip->updateFromKinetics(*kinetics);

  // update temperature dependent thermodynamic properties
  mythermo_->UpdateThermoProperties();

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

  return fromYAMLInput(pin->GetString("problem", input_key));
}

void Thermodynamics::Destroy() {
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (Thermodynamics::mythermo_ != nullptr) {
    delete Thermodynamics::mythermo_;
    Thermodynamics::mythermo_ = nullptr;
  }
}

size_t Thermodynamics::SpeciesIndex(std::string const& name) const {
  int index = kinetics_->thermo().speciesIndex(name);
  if (index < 0) {
    throw RuntimeError("Thermodynamics", "Species " + name + " not found");
  }
  return index;
}

std::string Thermodynamics::SpeciesName(size_t k) const {
  return kinetics_->thermo().speciesName(k);
}

Real Thermodynamics::GetTemp() const {
  return kinetics_->thermo().temperature();
}

Real Thermodynamics::GetPres() const { return kinetics_->thermo().pressure(); }

Real Thermodynamics::GetDensity() const {
  return kinetics_->thermo().density();
}

Real Thermodynamics::RovRd() const {
  auto& thermo = kinetics_->thermo();

  Real feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    feps += thermo.massFraction(n) * (inv_mu_ratio_[n] - 1.);
  }

  for (int n = 1 + NVAPOR; n < Size; ++n) {
    feps -= thermo.massFraction(n);
  }
  return feps;
}

void Thermodynamics::Extrapolate_inplace(Real dzORdlnp, std::string method,
                                         Real grav, Real userp) const {
  // RK4 integration
  if (grav == 0.) {  // hydrostatic
    _rk4_integrate_lnp(dzORdlnp, method, userp);
  } else {  // non-hydrostatic
    _rk4_integrate_z(dzORdlnp, method, grav, userp);
  }

  if (method != "reversible") {
    auto& thermo = kinetics_->thermo();
    std::vector<Real> xfrac(Size);
    thermo.getMoleFractions(xfrac.data());

    Real temp = thermo.temperature();
    Real pres = thermo.pressure();

    // a small number indicating saturation
    for (int j = 1 + NVAPOR; j < Size; ++j) {
      if (xfrac[j] > 0) xfrac[j] = 1.e-10;
    }
    thermo.setMoleFractions(xfrac.data());
    EquilibrateTP(temp, pres);
  }
}

Thermodynamics* Thermodynamics::mythermo_ = nullptr;

std::shared_ptr<Cantera::Condensation> get_kinetics_object(
    Thermodynamics const* pthermo) {
  return pthermo->kinetics_;
}
