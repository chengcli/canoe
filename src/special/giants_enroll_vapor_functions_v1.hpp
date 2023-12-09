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
#include <snap/thermodynamics/vapors/hydrogen_sulfide_vapors.hpp>
#include <snap/thermodynamics/vapors/methane_vapors.hpp>
#include <snap/thermodynamics/vapors/water_vapors.hpp>
#include <snap/thermodynamics/vapors/silicon_monoxide_vapors.hpp>

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

// silicon monoxide svp
void enroll_vapor_function_SiO(Thermodynamics::SVPFunc1Container& svp_func1,
                               std::vector<IndexSet> const& cloud_index_set) {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("SiO")) return;

  Application::Logger app("snap");
  app->Log("Enrolling SiO vapor pressures");

  int iSiO = pindex->GetVaporId("SiO");

  svp_func1[iSiO][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_SiO_Ferguson(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iSiO].size(); ++n) {
    svp_func1[iSiO][n] = NullSatVaporPres1;
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
  enroll_vapor_function_SiO(svp_func1_, cloud_index_set_);
  enroll_vapor_function_NH3(svp_func1_, cloud_index_set_);
  enroll_vapor_function_H2S(svp_func1_, cloud_index_set_);
  enroll_vapor_function_NH4SH(svp_func2_, cloud_reaction_map_);
}
