// C/C++ header
#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <index_map.hpp>
#include <variable.hpp>

// snap
#include "thermodynamics.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/ammonium_hydrosulfide_vapors.hpp"
#include "vapors/hydrogen_sulfide_vapors.hpp"
#include "vapors/methane_vapors.hpp"
#include "vapors/water_vapors.hpp"

Real NullSatVaporPres1(Variable const& qfrac, int i, int j) {
  Real g = 1.;

#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -qfrac.c[n];

  return qfrac.w[i] / g * qfrac.w[IPR];
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
      svp_func1_[iH2O][n] = [](Variable const& qfrac, int, int) {
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
      svp_func1_[iNH3][n] = [](Variable const& qfrac, int, int) {
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
      svp_func1_[iH2S][n] = [](Variable const& qfrac, int, int) {
        return sat_vapor_p_H2S_Antoine(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iH2S][n] = NullSatVaporPres1;
    }
  }
}

// methane svp
void Thermodynamics::enrollVaporFunctionCH4() {
  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("CH4")) return;

  Application::Logger app("snap");
  app->Log("Enrolling CH4 vapor pressures");

  int iCH4 = pindex->GetVaporId("CH4");
  for (int n = 0; n < cloud_index_set_[iCH4].size(); ++n) {
    int j = cloud_index_set_[iCH4][n];
    if (n == 0) {
      svp_func1_[iCH4][n] = [](Variable const& qfrac, int, int) {
        return sat_vapor_p_CH4_Antoine(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iCH4][n] = NullSatVaporPres1;
    }
  }
}

// ammonium hydrosulfide svp
void Thermodynamics::enrollVaporFunctionNH4SH() {
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
  svp_func2_[ij] = [](Variable const& qfrac, int, int, int) {
    return sat_vapor_p_NH4SH_Lewis(qfrac.w[IDN]);
  };

  cloud_reaction_map_[ij] = std::make_pair<ReactionIndx, ReactionStoi>(
      {iNH3, iH2S, iNH4SH}, {1, 1, 1});
}

void Thermodynamics::enrollVaporFunctionsGasGiants() {
  Application::Logger app("snap");
  app->Log("Enrolling vapor functions for gas giants");

  enrollVaporFunctionH2O();
  enrollVaporFunctionNH3();
  enrollVaporFunctionH2S();
  enrollVaporFunctionNH4SH();
}

void Thermodynamics::enrollVaporFunctionsIceGiants() {
  Application::Logger app("snap");
  app->Log("Enrolling vapor functions for ice giants");

  enrollVaporFunctionH2O();
  enrollVaporFunctionNH3();
  enrollVaporFunctionH2S();
  enrollVaporFunctionCH4();
  enrollVaporFunctionNH4SH();
}

void Thermodynamics::enrollVaporFunctionsEarth() { enrollVaporFunctionH2O(); }

void Thermodynamics::enrollVaporFunctionsMars() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsMars");
}

void Thermodynamics::enrollVaporFunctionsVenus() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsVenus");
}

void Thermodynamics::enrollVaporFunctionsTitan() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsTitan");
}

void Thermodynamics::enrollVaporFunctionsHotJupiter() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsHotJupiter");
}

void Thermodynamics::enrollVaporFunctionsLavaWorld() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsLavaWorld");
}
