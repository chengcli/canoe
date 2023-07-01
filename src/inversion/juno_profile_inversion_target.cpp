/** @file calculate_fit_target.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 21:07:59 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <iomanip>
#include <iostream>

// climath
#include <climath/interpolation.h>
#include <climath/linalg.h>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// inversion
#include "profile_inversion.hpp"

void JunoProfileInversion::CalculateFitTarget(Radiation const *prad, Real *val,
                                              int nvalue, int k, int j) const {
  std::stringstream msg;
  Application::Logger app("inversioon");

  app->Log("JunoProfileInversion::CalculateFitTarget");
  app->Log("model = " + std::to_string(j));

  if (nvalue != 2 * prad->GetNumBands()) {
    throw RuntimeError("CalculateFitTarget", "nvalue", prad->GetNumBands(),
                       nvalue);
  }

  // 11. log likelihood
  std::vector<Real> mus, tbs;

  for (int b = 0; b < prad->GetNumBands(); ++b) {
    auto pband = prad->GetBand(b);

    // emission angles;
    int ndir = pband->GetNumOutgoingRays();
    mus.resize(ndir);
    tbs.resize(ndir);

    for (int n = 0; n < ndir; ++n) mus[n] = pband->GetCosinePolarAngle(n);

    // brightness temperatures
    val[b * 2] = pband->btoa(0, k, j);

    // limb darkening
    for (int n = 0; n < ndir; ++n) tbs[n] = pband->btoa(n, k, j);

    Real tb45 = interp1(cos(45. / 180. * M_PI), tbs.data(), mus.data(), ndir);
    val[b * 2 + 1] = (tbs[0] - tb45) / tbs[0] * 100.;

    if (fit_differential_) {
      // brightness temperatures
      val[b * 2] -= pband->btoa(0, k, pmy_block_->js - 1);

      // limb darkening
      for (int n = 0; n < ndir; ++n)
        tbs[n] = pband->btoa(n, k, pmy_block_->js - 1);

      tb45 = interp1(cos(45. / 180. * M_PI), tbs.data(), mus.data(), ndir);
      val[b * 2 + 1] -= (tbs[0] - tb45) / tbs[0] * 100.;
    }
  }

  // app->Log("foward model results = ", val, nvalue);
}
