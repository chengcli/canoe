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
#include <snap/thermodynamics/molecules.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/vapors/ammonia_vapors.hpp>
#include <snap/thermodynamics/vapors/ammonium_hydrosulfide_vapors.hpp>
#include <snap/thermodynamics/vapors/hydrogen_sulfide_vapors.hpp>
#include <snap/thermodynamics/vapors/water_vapors.hpp>

// opacity
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/Giants/microwave/mwr_absorbers.hpp>

// utils
#include <utils/fileio.hpp>

// harp
#include <harp/radiation_band.hpp>
#include <opacity/absorber.hpp>

// outputs
#include <outputs/output_utils.hpp>

namespace gp = GiantPlanets;

// water svp
void enroll_vapor_function_H2O(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("H2O")) return;

  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  int iH2O = pindex->GetVaporId("H2O");

  svp_func1[iH2O][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_H2O_BriggsS(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iH2O].size(); ++n) {
    svp_func1[iH2O][n] = NullSatVaporPres1;
  }
}

// ammonia svp
void enroll_vapor_function_NH3(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("NH3")) return;

  Application::Logger app("snap");
  app->Log("Enrolling NH3 vapor pressures");

  int iNH3 = pindex->GetVaporId("NH3");

  svp_func1[iNH3][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_NH3_BriggsS(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iNH3].size(); ++n) {
    svp_func1[iNH3][n] = NullSatVaporPres1;
  }
}

// hydrogen sulfide svp
void enroll_vapor_function_H2S(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("H2S")) return;

  Application::Logger app("snap");
  app->Log("Enrolling H2S vapor pressures");

  int iH2S = pindex->GetVaporId("H2S");

  svp_func1[iH2S][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_H2S_Antoine(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iH2S].size(); ++n) {
    svp_func1[iH2S][n] = NullSatVaporPres1;
  }
}

// ammonium hydrosulfide svp
void enroll_vapor_function_NH4SH(
    Thermodynamics::SVPFunc2Container& svp_func2,
    std::map<IndexPair, ReactionInfo>& cloud_reaction_map) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("NH3") || !pindex->HasVapor("H2S") ||
      !pindex->HasCloud("NH4SH(s)"))
    return;

  Application::Logger app("snap");
  app->Log("Enrolling NH4SH vapor pressure");

  // ammonium hydrosulfide svp
  int iNH3 = pindex->GetVaporId("NH3");
  int iH2S = pindex->GetVaporId("H2S");
  int iNH4SH = pindex->GetCloudId("NH4SH(s)");

  auto ij = std::minmax(iNH3, iH2S);
  svp_func2[ij] = [](AirParcel const& qfrac, int, int, int) {
    return sat_vapor_p_NH4SH_Lewis(qfrac.w[IDN]);
  };

  cloud_reaction_map[ij] = std::make_pair<ReactionIndx, ReactionStoi>(
      {iNH3, iH2S, iNH4SH}, {1, 1, 1});
}

void Thermodynamics::enrollVaporFunctions() {
  Application::Logger app("snap");
  app->Log("Enrolling vapor functions");

  enroll_vapor_function_H2O(svp_func1_, cloud_index_set_);
  enroll_vapor_function_NH3(svp_func1_, cloud_index_set_);
  enroll_vapor_function_H2S(svp_func1_, cloud_index_set_);
  enroll_vapor_function_NH4SH(svp_func2_, cloud_reaction_map_);
}

extern Real xHe, xCH4;

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

// outputs
MetadataTable::MetadataTable() {
  Application::Logger app("outputs");
  app->Log("Initialize MetadataTable");

  table_ = {// short name, long name, units, grid location
            {"x1", "height at cell center", "m", "--C"},
            {"x1f", "height at cell boundary", "m", "--F"},
            {"x2", "distance at cell center", "m", "-C-"},
            {"x2f", "distance at cell boundary", "m", "-F-"},
            {"x3", "distance at cell center", "m", "C--"},
            {"x3f", "distance at cell boundary", "m", "F--"},
            {"rho", "density", "kg/m^3", "CCC"},
            {"press", "pressure", "pa", "CCC"},
            {"vel", "velocity", "m/s", "CCC"},
            {"vapor", "mass mixing ratio of vapor", "kg/kg", "CCC"},
            {"temp", "temperature", "K", "CCC"},
            {"theta", "potential temperature", "K", "CCC"},
            {"thetav", "virtual potential temperature", "K", "CCC"},
            {"mse", "moist static energy", "J/kg", "CCC"},
            {"rh1", "relative humidity 1", "1", "CCC"},
            {"rh2", "relative humidity 2", "1", "CCC"},
            {"eps", "turbulent dissipation", "w/kg", "CCC"},
            {"tke", "turbulent kinetic energy", "J/kg", "CCC"},
            {"mut", "dynamic turbulent viscosity", "kg/(m.s)", "CCC"},
            {"CH1tau", "optical thickness", "1", "CCC"},
            {"CH2tau", "optical thickness", "1", "CCC"},
            {"CH3tau", "optical thickness", "1", "CCC"},
            {"CH4tau", "optical thickness", "1", "CCC"},
            {"CH5tau", "optical thickness", "1", "CCC"},
            {"CH6tau", "optical thickness", "1", "CCC"},
            {"radiance", "top-of-atmosphere radiance", "K", "RCC"}};
}
