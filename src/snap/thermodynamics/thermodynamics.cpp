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

  for (int n = 1 + NVAPOR; n < Size; ++n) {
    cv_ratio_[n] = gammad_ * cp_ratio_[n];
  }
}

Thermodynamics* Thermodynamics::fromYAMLInput(std::string const& fname) {
  mythermo_ = new Thermodynamics();

  auto thermo = Cantera::newThermo(fname, "gas");

  if (thermo->nSpecies() != 1 + NVAPOR + NCLOUD) {
    throw RuntimeError("Thermodynamics",
                       "Number of species does not match the input file");
  }

  auto& kinetics = mythermo_->kinetics_;
  kinetics = std::static_pointer_cast<Cantera::Condensation>(
      Cantera::newKinetics({thermo}, fname));
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

  return fromYAMLInput(pin->GetString("thermodynamics", input_key));
}

void Thermodynamics::Destroy() {
  std::unique_lock<std::mutex> lock(thermo_mutex);

  if (Thermodynamics::mythermo_ != nullptr) {
    delete Thermodynamics::mythermo_;
    Thermodynamics::mythermo_ = nullptr;
  }
}

Thermodynamics* Thermodynamics::mythermo_ = nullptr;
