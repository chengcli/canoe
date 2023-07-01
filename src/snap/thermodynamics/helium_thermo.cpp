/** @file helium.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday Jul 13, 2021 15:26:31 PDT
 * @bug No known bugs.
 */

#include <sstream>
#include <stdexcept>

#include "molecules.hpp"

double Helium::nist_shomate_[7] = {20.78603,     4.850638E-10, -1.582916E-10,
                                   1.525102E-11, 3.196347E-11, -6.197341,
                                   151.3064};

double Helium::cp_nist(double T) {
  std::stringstream msg;
  T = std::min(std::max(298., T), 600.);
  // if (T < 298. || T > 600.) {
  //   msg << "ERROR: Temperature out of range in Helium::cp_nist" << std::endl;
  //   throw std::runtime_error(msg.str().c_str());
  // }

  double *pdata = nist_shomate_;
  double result;
  T /= 1.E3;
  result = pdata[0] + pdata[1] * T + pdata[2] * T * T + pdata[3] * T * T * T +
           pdata[4] / (T * T);
  return result;
}

Helium aHe;
