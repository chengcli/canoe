// C/C++
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <index_map.hpp>

// thermodynamics
#include <thermodynamics/vapors/ammonia_vapors.hpp>
#include <thermodynamics/vapors/ammonium_hydrosulfide_vapors.hpp>
#include <thermodynamics/vapors/hydrogen_sulfide_vapors.hpp>
#include <thermodynamics/vapors/methane_vapors.hpp>
#include <thermodynamics/vapors/water_vapors.hpp>

// opacity
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/Giants/microwave/mwr_absorbers.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/parameter_map.hpp>

// harp
#include <harp/absorber.hpp>
#include <harp/radiation_band.hpp>

namespace gp = GiantPlanets;

// thermodynamics
Real NullSatVaporPres1(AirParcel const& qfrac, int i, int j) {
  Real g = 1.;

#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -qfrac.c[n];

  return qfrac.w[i] / g * qfrac.w[IPR];
}

void Thermodynamics::enrollVaporFunctions() {
  Application::Logger app("snap");
  app->Log("Enrolling vapor functions for juno mwr");

  enrollVaporFunctionH2O();
  enrollVaporFunctionNH3();
  enrollVaporFunctionH2S();
  enrollVaporFunctionNH4SH();
}

// water svp
void Thermodynamics::enrollVaporFunctionH2O() {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("H2O")) return;

  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  int iH2O = pindex->GetVaporId("H2O");
  for (int n = 0; n < cloud_index_set_[iH2O].size(); ++n) {
    int j = cloud_index_set_[iH2O][n];
    if (n == 0) {
      svp_func1_[iH2O][n] = [](AirParcel const& qfrac, int, int) {
        return sat_vapor_p_H2O_BriggsS(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iH2O][n] = NullSatVaporPres1;
    }
  }
}

// ammonia svp
void Thermodynamics::enrollVaporFunctionNH3() {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("NH3")) return;

  Application::Logger app("snap");
  app->Log("Enrolling NH3 vapor pressures");

  int iNH3 = pindex->GetVaporId("NH3");
  for (int n = 0; n < cloud_index_set_[iNH3].size(); ++n) {
    int j = cloud_index_set_[iNH3][n];
    if (n == 0) {
      svp_func1_[iNH3][n] = [](AirParcel const& qfrac, int, int) {
        return sat_vapor_p_NH3_BriggsS(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iNH3][n] = NullSatVaporPres1;
    }
  }
}

// hydrogen sulfide svp
void Thermodynamics::enrollVaporFunctionH2S() {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("H2S")) return;

  Application::Logger app("snap");
  app->Log("Enrolling H2S vapor pressures");

  int iH2S = pindex->GetVaporId("H2S");
  for (int n = 0; n < cloud_index_set_[iH2S].size(); ++n) {
    int j = cloud_index_set_[iH2S][n];
    if (n == 0) {
      svp_func1_[iH2S][n] = [](AirParcel const& qfrac, int, int) {
        return sat_vapor_p_H2S_Antoine(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iH2S][n] = NullSatVaporPres1;
    }
  }
}

// hydrogen heat capacity
Real Thermodynamics::GetGammad(AirParcel const& qfrac) const {
  Real T = qfrac.w[IDN], cp_h2, cp_he, cp_ch4;
  if (T < 300.) {
    cp_h2 = Hydrogen::cp_norm(T);
  } else {
    cp_h2 = Hydrogen::cp_nist(T);
  }
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

  Real cp_real = (1. - xHe - xCH4) * cp_h2 + xHe * cp_he + xCH4 * cp_ch4;
  return cp_real / (cp_real - Constants::Rgas);
}

// radiation
void RadiationBand::addAbsorberRadio(YAML::Node& node) {
  auto name = node["name"].as<std::string>();
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "NH3") {
    auto ab = std::make_unique<gp::MwrAbsorberNH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2O") {
    auto ab = std::make_unique<gp::MwrAbsorberH2O>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2S") {
    auto ab = std::make_unique<gp::MwrAbsorberH2S>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "PH3") {
    auto ab = std::make_unique<gp::MwrAbsorberPH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CIA") {
    auto ab = std::make_unique<gp::MwrAbsorberCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "Electron") {
    auto ab = std::make_unique<gp::MwrAbsorberElectron>(species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberRadio", "Absorber " + name);
  }
}
