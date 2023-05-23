#ifndef WATER_VAPORS_HPP_
#define WATER_VAPORS_HPP_

#include <cmath>

inline double svph2o(double t, double p, double beta, double gamma) {
  return p * exp((1. - 1. / t) * beta - gamma * log(t));
}

inline double sat_vapor_p_H2O_UMich(double T) {
  double x = -2445.5646 / T + 8.2312 * log10(T) - 0.01677006 * T +
             1.20514e-05 * T * T - 6.757169;
  return pow(10.0, x) * 1.E5 / 760.;
}

// (273, 373)
inline double sat_vapor_p_H2O_Antoine(double T) {
  double result;

  if (T < 303.)
    result = pow(10., 5.40221 - (1838.675 / (T - 31.737)));
  else if (T < 333.)
    result = pow(10., 5.20389 - (1733.926 / (T - 39.485)));
  else if (T < 363.)
    result = pow(10., 5.0768 - (1659.793 / (T - 45.854)));
  else
    result = pow(10., 5.08354 - (1663.125 / (T - 45.622)));

  return 1.E5 * result;
}

// (100, 373.16)
inline double sat_vapor_p_H2O_Hubner(double T) {
  double A, B, C, D;
  A = 4.07023;
  B = -2484.98;
  C = 3.56654;
  D = -0.00320981;
  double x = A + B / T + C * log10(T) + D * T;  // best fit in [100K; 273.16K]
  return pow(10.0, x);
}

// (0, 273.16)
inline double sat_vapor_p_H2O_Fray(double T) {
  double a[7], x, pt, Tt, tr;
  pt = 6.11657e-03;
  Tt = 273.16;
  tr = T / Tt;
  x = 0.;
  a[0] = 20.9969665107897;
  a[1] = 3.72437478271362;
  a[2] = -13.9205483215524;
  a[3] = 29.6988765013566;
  a[4] = -40.1972392635944;
  a[5] = 29.7880481050215;
  a[6] = -9.13050963547721;

  for (int i = 0; i < 7; i++) x = x + a[i] * pow(tr, i);
  return 1.E5 * pt * exp(1.5 * log(tr) + (1 - 1 / tr) * x);
}

inline double sat_vapor_p_H2O_BriggsS(double T) {
  double a[6], x;
  if (T < 273.16) {
    a[1] = -5631.1206;
    a[2] = -8.363602;
    a[3] = 8.2312;
    a[4] = -3.861449e-02;
    a[5] = 2.77494e-05;
  } else {
    a[1] = -2313.0338;
    a[2] = -164.03307;
    a[3] = 38.053682;
    a[4] = -1.3844344e-01;
    a[5] = 7.4465367e-05;
  }
  x = a[1] / T + a[2] + a[3] * log(T) + a[4] * T + a[5] * pow(T, 2);
  return exp(x) / 10.;
}

inline double sat_vapor_p_H2O_Bolton(double T) {
  double result;
  result = 612.2 * exp(17.67 * (T - 273.15) / (T - 29.65));
  return result;
}

inline double sat_vapor_p_H2O_Smithsonian(double T) {
  double result;
  result = 100. * exp(23.33086 - 6111.72784 / T + 0.15215 * log(T));
  return result;
}

inline double sat_vapor_p_H2O_Ideal(double T) {
  double betal = 24.88, gammal = 5.06, betas = 22.98, gammas = 0.52,
         tr = 273.16, pr = 611.7;
  return T > tr ? svph2o(T / tr, pr, betal, gammal)
                : svph2o(T / tr, pr, betas, gammas);
}

#endif
