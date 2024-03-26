#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_

// C/C++
#include <cmath>

// best fit in [154.26, 195.89] K
inline double sat_vapor_p_CO2_Antoine(double T) {
  double A = 6.81228;
  double B = 1301.679;
  double C = -34.94;

  double result = pow(10, A - B / (T + C));
  return 1.E5 * result;

  // best fit in [154.2, 204] K

  // double A = 9.8106;
  // double B = 1347.79;
  // double C = 273;

  // double result = pow(10, A - B / ((T-273.2) + C));
  // return 133.3 * result;
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_CARBON_DIOXIDE_VAPORS_HPP_
