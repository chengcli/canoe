#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_SULFUR_DIOXIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_SULFUR_DIOXIDE_VAPORS_HPP_

// C/C++
#include <cmath>

// best fit in [177.7, 263] K
inline double sat_vapor_p_SO2_Antoine(double T) {
  double A = 3.48586;
  double B = 668.225;
  double C = -72.252;
  double result = pow(10, A - B / (T + C));
  return 1.E5 * result;
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_SULFUR_DIOXIDE_VAPORS_HPP_
