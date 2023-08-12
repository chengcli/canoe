// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// snap
#include "thermodynamics.hpp"

Real NullSatVaporPres1(AirParcel const& air, int i, int j) {
  Real g = 1.;

#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -air.c[n];

  return air.w[i] / g * air.w[IPR];
}

// all vapor functions
void __attribute__((weak)) Thermodynamics::enrollVaporFunctions() {}

// water svp
void __attribute__((weak)) Thermodynamics::enrollVaporFunctionH2O() {}

// ammonia svp
void __attribute__((weak)) Thermodynamics::enrollVaporFunctionNH3() {}

// hydrogen sulfide svp
void __attribute__((weak)) Thermodynamics::enrollVaporFunctionH2S() {}

// methane svp
void __attribute__((weak)) Thermodynamics::enrollVaporFunctionCH4() {}

// ammonium hydrosulfide svp
void __attribute__((weak)) Thermodynamics::enrollVaporFunctionNH4SH() {}
