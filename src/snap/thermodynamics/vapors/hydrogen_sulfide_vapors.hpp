#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_HYDROGEN_SULFIDE_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_HYDROGEN_SULFIDE_VAPORS_HPP_

// C/C++
#include <cmath>

double Pcgs_from_torr(double P) { return (1013250.0 / 760.0) * P; }

double sat_vapor_p_H2S_solid_UMich(double T) {
  //  Vapor of H2S over H2S liquid, T > 188 K, over H2S ice for T < 188 K
  double x;
  double log10T = log10(T);
  if (T > 212.75)
    x = 7.7547 - 976.0 / T - 0.12058 * log10T;
  else if (T > 188.0)
    x = 15.4859 - 1264.3574 / T - 2.86206 * log10T;
  else
    x = 57.19 - 2461.84 / T - 18.443 * log10T;

  return Pcgs_from_torr(pow(10.0, x)) * 0.1;
}

//   best fit in [138.8K; 349.5K] (webbook.nist.gov)
double sat_vapor_p_H2S_Antoine(double T) {
  double result;
  if (T < 212.8)
    result = pow(10., 4.43681 - (829.439 / (T - 25.412)));
  else
    result = pow(10., 4.52887 - (958.587 / (T - 0.539)));

  return 1.E5 * result;
}

// best fit in [120K; 190K]
double sat_vapor_p_H2S_solid_Hubner(double T) {
  double A, B, C, D;
  double log10T = log10(T);
  A = 6.96156;
  B = -903.815;
  C = 0.258812;
  D = 0.00873804;

  double x = A + B / T + C * log10(T) + D * T;
  return 1.E5 * pow(10.0, x);
}

// best fit in [110K; 187.57K]
double sat_vapor_p_H2S_solid_Fray(double T) {
  double a[7], x;
  x = 0.;
  if (T > 127) {
    a[0] = 8.933;
    a[1] = -7.260e+02;
    a[2] = -3.504e+05;
    a[3] = 2.724e+07;
    a[4] = -8.582e+08;
    a[5] = 0.;
    a[6] = 0.;
  } else {
    a[0] = 1.298e+01;
    a[1] = -2.707e+03;
    a[2] = 0.;
    a[3] = 0.;
    a[4] = 0.;
    a[5] = 0.;
    a[6] = 0.;
  }

  for (int i = 1; i < 7; i++) {
    x = x + a[i] / pow(T, i);
  }
  return 1.E5 * exp(x + a[0]);
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_HYDROGEN_SULFIDE_VAPORS_HPP_
