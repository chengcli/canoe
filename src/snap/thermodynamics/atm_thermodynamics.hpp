#pragma once

#include "thermodynamics.hpp"

template <typename T>
T potential_temperature(Thermodynamics *pthermo, T w, float p0) {
  return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
}

template <typename T>
T equivalent_potential_temperature(Thermodynamics *pthermo, T w, float p0) {
  return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
}

template <typename A, typename B>
A moist_static_energy(Thermodynamics *pthermo, A w, B gz) {
  return pthermo->MoistEnthalpy(w) + gz;
}

template <typename T>
T relative_humidity(Thermodynamics *pthermo, T w, int ivapor) {
  return pthermo->GetDensity(iH2O) / w[IDN];
}
