// C/C++ header
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
#include "vapors/water_vapors.hpp"

Real NullSatVaporPres(Variable const& qfrac, int i, int j) {
  Real g = 1.;

#pragma omp reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g -= qfrac.c[n];

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
      svp_func_[iH2O][n] = [](Variable const& qfrac, int, int) {
        return sat_vapor_p_H2O_BriggsS(qfrac.w[IDN]);
      };
    } else {
      svp_func_[iH2O][n] = NullSatVaporPres;
    }
  }

  // ammonia svp:
  int iNH3 = pindex->GetVaporId("NH3");
  for (int n = 0; n < cloud_index_set_[iNH3].size(); ++n) {
    int j = cloud_index_set_[iNH3][n];
    if (n == 0) {
      svp_func_[iNH3][n] = [](Variable const& qfrac, int, int) {
        return sat_vapor_p_NH3_BriggsS(qfrac.w[IDN]);
      };
    } else {
      svp_func_[iNH3][n] = NullSatVaporPres;
    }
  }
}

void Thermodynamics::enrollVaporFunctionsEarth() {
  throw NotImplementedError("Thermodynamics::enrollVaporFunctionsEarth");
}
