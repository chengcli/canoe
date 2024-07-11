// application
#include <application/application.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <index_map.hpp>

// thermodynamics
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/vapors/ammonia_vapors.hpp>
#include <snap/thermodynamics/vapors/ammonium_hydrosulfide_vapors.hpp>
#include <snap/thermodynamics/vapors/carbon_dioxide_vapors.hpp>
#include <snap/thermodynamics/vapors/hydrogen_sulfide_vapors.hpp>
#include <snap/thermodynamics/vapors/methane_vapors.hpp>
#include <snap/thermodynamics/vapors/sulfur_dioxide_vapors.hpp>
#include <snap/thermodynamics/vapors/water_vapors.hpp>

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

// sulfur dioxide svp
void enroll_vapor_function_SO2(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("SO2")) return;

  Application::Logger app("snap");
  app->Log("Enrolling SO2 vapor pressures");

  int iSO2 = pindex->GetVaporId("SO2");

  svp_func1[iSO2][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_SO2_Antoine(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iSO2].size(); ++n) {
    svp_func1[iSO2][n] = NullSatVaporPres1;
  }
}

// carbon dioxide svp
void enroll_vapor_function_CO2(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("CO2")) return;

  Application::Logger app("snap");
  app->Log("Enrolling CO2 vapor pressures");

  int iCO2 = pindex->GetVaporId("CO2");

  svp_func1[iCO2][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_CO2_Antoine(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iCO2].size(); ++n) {
    svp_func1[iCO2][n] = NullSatVaporPres1;
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
  enroll_vapor_function_SO2(svp_func1_, cloud_index_set_);
  enroll_vapor_function_CO2(svp_func1_, cloud_index_set_);
  enroll_vapor_function_NH4SH(svp_func2_, cloud_reaction_map_);
}
