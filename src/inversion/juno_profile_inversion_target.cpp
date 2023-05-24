/** @file calculate_fit_target.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 21:07:59 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <climath/interpolation.h>
#include <climath/linalg.h>

#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <debugger/debugger.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>
#include <iomanip>
#include <iostream>

#include "profile_inversion.hpp"

extern std::unique_ptr<Debugger> pdebug;

void JunoProfileInversion::CalculateFitTarget(Radiation const *prad, Real *val,
                                              int nvalue, int k, int j) const {
  std::stringstream msg;
  pdebug->Call("JunoProfileInversion::CalculateFitTarget");
  pdebug->Message("model", j);

  // 11. log likelihood
  std::vector<Real> mus, tbs;

  int b = 0, bid = 0;
  for (auto p : prad->bands) {
    // emission angles;
    int ndir = p->getNumOutgoingRays();
    mus.resize(ndir);
    tbs.resize(ndir);

    for (int n = 0; n < ndir; ++n) mus[n] = p->getCosinePolarAngle(n);

    // brightness temperatures
    val[b * 2] = prad->radiance(bid, k, j);

    // limb darkening
    for (int n = 0; n < ndir; ++n) tbs[n] = prad->radiance(bid + b, k, j);

    Real tb45 = interp1(cos(45. / 180. * M_PI), tbs.data(), mus.data(), ndir);
    val[b * 2 + 1] = (tbs[0] - tb45) / tbs[0] * 100.;

    if (fit_differential_) {
      // brightness temperatures
      val[b * 2] -= prad->radiance(bid, k, pmy_block_->js - 1);

      // limb darkening
      for (int n = 0; n < ndir; ++n)
        tbs[n] = prad->radiance(bid, k, pmy_block_->js - 1);

      tb45 = interp1(cos(45. / 180. * M_PI), tbs.data(), mus.data(), ndir);
      val[b * 2 + 1] -= (tbs[0] - tb45) / tbs[0] * 100.;
    }

    bid += ndir;
    b++;
  }

  pdebug->Message("foward model results", val, nvalue);

  pdebug->Leave();
}
