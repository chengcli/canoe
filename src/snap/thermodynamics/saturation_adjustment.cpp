// C/C++
#include <cstdlib>
#include <cstring>
#include <sstream>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/hydro/hydro.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::SaturationAdjustment(AirColumn& air_column) const {
  // return if there's no vapor
  if (NVAPOR == 0) return;

  for (auto& air : air_column) {
    EquilibrateUV(&air);
  }
}
