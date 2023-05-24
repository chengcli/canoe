#ifndef SRC_SNAP_THERMODYNAMICS_VAPORS_SODIUM_VAPORS_HPP_
#define SRC_SNAP_THERMODYNAMICS_VAPORS_SODIUM_VAPORS_HPP_

// TODO(cli) check this equation
inline double sat_vapor_p_Na_H2S_Visscher(double T, double pH2S) {
  double log10p = 8.55 - 13889. / T - 0.5 * log10(pH2S / 1E5);
  return 1.E5 * pow(10., log10p);
}

#endif  // SRC_SNAP_THERMODYNAMICS_VAPORS_SODIUM_VAPORS_HPP_
