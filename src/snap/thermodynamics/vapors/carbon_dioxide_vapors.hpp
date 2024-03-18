#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_

// C/C++
#include <cmath>

// best fit in [154.26, 195.89] K
inline double sat_vapor_p_CO2_Antoine(double T) {
  double A = 6.81228;
  double B = 1301.679;
  double C = -3.494;
  double result = pow(10, A - B / (T + C));
  return 1.E5 * result;
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_
