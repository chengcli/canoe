// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// snap
#include "thermodynamics.hpp"

Real NullSatVaporPres1(AirParcel const& air, int i, int j) {
  Real g = 1.;

<<<<<<< HEAD
#pragma omp simd reduction(+ : g)
  for (int n = 0; n < NCLOUD; ++n) g += -air.c[n];
=======
#pragma omp simd reduction(- : g)
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
>>>>>>> 027e882 (intel)

  return air.w[i] / g * air.w[IPR];
}

// enroll null vapor functions
void __attribute__((weak)) Thermodynamics::enrollVaporFunctions() {
  for (int i = 0; i <= NVAPOR; ++i)
    for (int j = 0; j < NPHASE - 1; ++j) {
      svp_func1_[i][j] = NullSatVaporPres1;
    }
}
