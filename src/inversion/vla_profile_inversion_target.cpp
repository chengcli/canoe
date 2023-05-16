/** @file calculate_fit_target.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 21:07:59 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>
#include <iomanip>

// Athena++ headers
#include <athena.hpp>
#include <mesh/mesh.hpp>

// climath headers
extern "C" {
  #include <linalg.h>
  #include <interpolation.h>
}

// harp2 headers
#include "../radiation/radiation.hpp"
#include "../radiation/radiation_band.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/block_index.hpp"
#include "profile_inversion.hpp"

extern std::unique_ptr<Debugger> pdebug;

void VLAProfileInversion::CalculateFitTarget(Radiation const *prad,
  Real *val, int nvalue, int k, int j) const
{
  std::stringstream msg;
  pdebug->Call("VLAProfileInversion::CalculateFitTarget");
  pdebug->Message("model", j);

  // 11. log likelihood
  std::vector<Real> mus, tbs;

  int b = 0, bid = 0;
  for (auto p : prad->bands) {
    // emission angles;
    int ndir = p->getNumOutgoingRays();
    mus.resize(ndir);
    tbs.resize(ndir);

    for (int n = 0; n < ndir; ++n)
      mus[n] = p->getCosinePolarAngle(n);

    // brightness temperature
    val[b] = prad->radiance(bid,k,j);

    if (fit_differential_) {
      // brightness temperature differential
      val[b] -= prad->radiance(bid,k,pblock_->js-1);
    }

    bid += ndir;
    b++;
    if (b >= nvalue) break;
  }

  pdebug->Message("foward model results", val, nvalue);

  pdebug->Leave();
}
