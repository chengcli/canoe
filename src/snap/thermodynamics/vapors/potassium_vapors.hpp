/** @file potassium_vapors.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Dec 04, 2021 22:40:24 EST
 * @bug No known bugs.
 */

#ifndef POTASSIUM_VAPORS_HPP
#define POTASSIUM_VAPORS_HPP

inline double sat_vapor_p_KCl_Lodders(double T) {
  double logp = 7.611 - 11382. / T;
  return 1.E5 * exp(logp);
}

#endif /* end of include guard POTASSIUM_VAPORS_HPP */
