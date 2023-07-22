// C/C++ header
#include <algorithm>
#include <fstream>
#include <iostream>

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
#include "vapors/water_vapors.hpp"

Real NullSatVaporPres1(Variable const& qfrac, int i, int j) {
  Real g = 1.;

#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -qfrac.c[n];

  return qfrac.w[i] / g * qfrac.w[IPR];
}

void Thermodynamics::enrollVaporFunctions(ParameterInput* pin) {
  if (NVAPOR == 0 && NCLOUD == 0) return;

  if (strcmp(PLANET, "Earth") == 0) {
    enrollVaporFunctionsEarth();
  } else if (strcmp(PLANET, "Jupiter") == 0) {
    enrollVaporFunctionsGiants();
  } else {
    throw RuntimeError("Thermodynamics",
                       "Unknown planet: " + std::string(PLANET));
  }
}

void Thermodynamics::enrollVaporFunctionsGiants() {
  enrollVaporFunctionsJupiterJuno();
}

void Thermodynamics::enrollVaporFunctionsJupiterJuno() {
  auto pindex = IndexMap::GetInstance();

  // water svp:
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

  // ammonia svp:
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

  // hydrogen sulfide svp:
  int iH2S = pindex->GetVaporId("H2S");
  for (int n = 0; n < cloud_index_set_[iH2S].size(); ++n) {
    int j = cloud_index_set_[iH2S][n];
    if (n == 0) {
      svp_func1_[iH2S][n] = [](Variable const& qfrac, int, int) {
        return sat_vapor_p_H2S_BriggsS(qfrac.w[IDN]);
      };
    } else {
      svp_func1_[iH2S][n] = NullSatVaporPres1;
    }
  }

  // ammonium hydrosulfide svp
  int iNH4SH = pindex->GetCloudId("NH4SH(s)");
  auto ij = std::minmax(iNH3, iH2S);
  svp_func2_[ij] = [](Variable const& qfrac, IndexPair _ij) {
    return sat_vapor_p_NH4SH(qfrac.w[IDN]);
  };
  cloud_reaction_map_[ij] = std::make_pair({iNH3, iH2S, iNH4SH}, {1, 1, 1});
}

void Thermodynamics::enrollVaporFunctionsEarth() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsEarth");
}
