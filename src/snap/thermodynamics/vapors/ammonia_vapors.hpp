#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIA_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIA_VAPORS_HPP_

#include <cmath>

inline double svpnh3(double t, double p, double beta, double gamma) {
  return p * exp((1. - 1. / t) * beta - gamma * log(t));
}

inline double sat_vapor_p_NH3_UMich(double T) {
  double x = -1790.00 / T - 1.81630 * log10(T) + 14.97593;
  return pow(10.0, x) * 1.E5 / 760.;
}

// (164, 371.5)
inline double sat_vapor_p_NH3_Antoine(double T) {
  double result;
  if (T < 239.6)
    result = pow(10., 3.18757 - (506.713 / (T - 80.78)));
  else
    result = pow(10., 4.86886 - (1113.928 / (T - 10.409)));
  return 1.E5 * result;
}

// (130, 200)
inline double sat_vapor_p_NH3_Hubner(double T) {
  double A = 24.3037, B = -1766.28, C = -5.64472, D = 0.00740241;

  double x = A + B / T + C * log10(T) + D * T;
  return pow(10.0, x);
}

inline double sat_vapor_p_NH3_BriggsS(double T) {
  double a[6], x;
  if (T < 195) {
    a[1] = -4122.;
    a[2] = 41.67871;
    a[3] = -1.81630;
    a[4] = 0.;
    a[5] = 0.;
  } else {
    a[1] = -4409.3512;
    a[2] = 76.864252;
    a[3] = -8.4598340;
    a[4] = 5.51029e-03;
    a[5] = 6.80463e-06;
  }
  x = a[1] / T + a[2] + a[3] * log(T) + a[4] * T + a[5] * pow(T, 2);
  return exp(x) / 10.;
}

// (15, 195.4)
inline double sat_vapor_p_NH3_Fray(double T) {
  double a[7], x = 0;
  a[0] = 1.596e+01;
  a[1] = -3.537e+03;
  a[2] = -3.310e+04;
  a[3] = 1.742e+06;
  a[4] = -2.995e+07;
  a[5] = 0.;
  a[6] = 0.;

  for (int i = 1; i < 7; i++)
    x = x + a[i] / pow(T, i);  // best fit in [15K; 195.41K]

  return 1.E5 * exp(x + a[0]);
}

inline double sat_vapor_p_NH3_Ideal(double T) {
  double betal = 20.08, gammal = 5.62, betas = 20.64, gammas = 1.43, tr = 195.4,
         pr = 6060.;

  return T > tr ? svpnh3(T / tr, pr, betal, gammal)
                : svpnh3(T / tr, pr, betas, gammas);
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_AMMONIA_VAPORS_HPP_
