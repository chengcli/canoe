#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_METHANE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_METHANE_VAPORS_HPP_

// C/C++
#include <cmath>

inline double svpch4(double t, double p, double beta, double gamma) {
  return p * exp((1. - 1. / t) * beta - gamma * log(t));
}

inline double sat_vapor_p_CH4_Ideal(double T) {
  double betal = 10.15, gammal = 2.1;
  double betas = 10.41, gammas = 0.9;
  double tr = 90.67, pr = 11690.;

  return T > tr ? svpch4(T / tr, pr, betal, gammal)
                : svpch4(T / tr, pr, betas, gammas);
}

// best fit in [90.99, 189.99] K
inline double sat_vapor_p_CH4_Antoine(double T) {
  double A = 3.9895;
  double B = 443.028;
  double C = -0.49;
  double result = pow(10, A - B / (T + C));
  return 1.E5 * result;
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_METHANE_VAPORS_HPP_
