// canoe
#include <configure.h>

#include <air_parcel.hpp>

// snap
#include "thermodynamics.hpp"

Real NullSatVaporPres1(AirParcel const& air, int i, int j) {
  Real g = 1.;

#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -air.c[n];

  return air.w[i] / g * air.w[IPR];
}

// enroll null vapor functions
void __attribute__((weak)) Thermodynamics::enrollVaporFunctions() {
  for (int i = 0; i <= NVAPOR; ++i)
    for (int j = 0; j < NPHASE - 1; ++j) {
      svp_func1_[i][j] = NullSatVaporPres1;
    }
}
