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

// inversion
#include <application/application.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// inversion
#include "profile_inversion.hpp"

void VLAProfileInversion::CalculateFitTarget(Radiation const *prad, Real *val,
                                             int nvalue, int k, int j) const {
  std::stringstream msg;
  Application::Logger app("inversion");

  app->Log("VLAProfileInversion::CalculateFitTarget");
  app->Log("model = " + std::to_string(j));

  // 11. log likelihood
  std::vector<Real> mus, tbs;

  int b = 0, bid = 0;
  for (auto p : prad->bands) {
    // emission angles;
    int ndir = p->GetNumOutgoingRays();
    mus.resize(ndir);
    tbs.resize(ndir);

    for (int n = 0; n < ndir; ++n) mus[n] = p->getCosinePolarAngle(n);

    // brightness temperature
    val[b] = prad->radiance(bid, k, j);

    if (fit_differential_) {
      // brightness temperature differential
      val[b] -= prad->radiance(bid, k, pmy_block_->js - 1);
    }

    bid += ndir;
    b++;
    if (b >= nvalue) break;
  }

  // app->Log("foward model results", val, nvalue);
}
