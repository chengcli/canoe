#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIUM_HYDROSULFIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIUM_HYDROSULFIDE_VAPORS_HPP_

// C/C++
#include <cmath>

const double Pcgs_of_atm = 1013250.0;  // atmospheres to dynes/cm**2

inline double sat_vapor_p_NH4SH_UMich(double T) {
  double const GOLB2 = (14.83 - (4715.0 / T));
  return (pow(10.0, GOLB2)) * Pcgs_of_atm * Pcgs_of_atm;
}

inline double sat_vapor_p_NH4SH_Lewis(double T) {
  return pow(10., 14.82 - 4705. / T) * 101325. * 101325.;
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIUM_HYDROSULFIDE_VAPORS_HPP_
