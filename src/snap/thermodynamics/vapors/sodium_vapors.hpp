/** @file sodium_vapors.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Dec 04, 2021 22:43:54 EST
 * @bug No known bugs.
 */

#ifndef SODIUM_VAPORS_HPP
#define SODIUM_VAPORS_HPP

// TODO check this equation
inline double sat_vapor_p_Na_H2S_Visscher(double T, double pH2S) {
  double log10p = 8.55 - 13889. / T - 0.5 * log10(pH2S / 1E5);
  return 1.E5 * pow(10., log10p);
}

#endif /* end of include guard SODIUM_VAPORS_HPP */
