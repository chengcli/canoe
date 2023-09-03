/** @file methane.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday Jul 13, 2021 15:26:31 PDT
 * @bug No known bugs.
 */

#include <sstream>
#include <stdexcept>

#include "molecules.hpp"

double Methane::nist_shomate1_[7] = {-0.703029, 108.4773,  -42.52157, 5.862788,
                                     0.678565,  -76.84376, 158.7163};

double Methane::nist_shomate2_[7] = {85.81217,  11.26467,  -2.114146, 0.138190,
                                     -26.42221, -153.5327, 224.4143};

double Methane::cp_nist(double T) {
  double *pdata;
  std::stringstream msg;
  T = std::min(std::max(298., T), 6000.);
  // if (T < 298. || T > 6000.) {
  //   msg << "ERROR: Temperature out of range in Methane::cp_nist" <<
  //   std::endl; throw std::runtime_error(msg.str().c_str());
  // }

  if (T < 1300.) {
    pdata = nist_shomate1_;
  } else {
    pdata = nist_shomate2_;
  }

  double result;
  T /= 1.E3;
  result = pdata[0] + pdata[1] * T + pdata[2] * T * T + pdata[3] * T * T * T +
           pdata[4] / (T * T);
  return result;
}

Methane aCH4;
