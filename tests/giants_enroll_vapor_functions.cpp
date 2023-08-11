// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <index_map.hpp>

// thermodynamics
#include <snap/thermodynamics/thermodynamcis.hpp>
#include <thermodynamics/vapors/ammonia_vapors.hpp>
#include <thermodynamics/vapors/ammonium_hydrosulfide_vapors.hpp>
#include <thermodynamics/vapors/hydrogen_sulfide_vapors.hpp>
#include <thermodynamics/vapors/methane_vapors.hpp>
#include <thermodynamics/vapors/water_vapors.hpp>

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
