#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_SILICON_MONOXIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_SILICON_MONOXIDE_VAPORS_HPP_

// C/C++
#include <cmath>

inline double sat_vapor_p_SiO_Ferguson(double T) {
  double log10p = 13.25 - 17900. / T;
  return pow(10., log10p);
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_SILICON_MONOXIDE_VAPORS_HPP_
