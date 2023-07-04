// C/C++
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// canoe
#include <configure.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/vectorize.hpp>

// thermodynamics
#include "thermodynamics_helper.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/water_vapors.hpp"

void __attribute__((weak)) update_gamma(Real *gamma, Real const q[]) {}

/*Real qhat_eps(Real const q[], Real const eps[])
{
  Real feps = 1.;
  for (int n = 1; n < NMASS; ++n)
    feps += q[n]*(eps[n] - 1.);
  return q_gas(q)/feps;
}*/
