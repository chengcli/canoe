/** @file pi_log_posterior_probability.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Nov 20, 2021 09:00:13 EST
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>
#include <memory>
#include <sstream>

// athena
#include <athena/mesh/mesh.hpp>

// harp
#include <harp/radiation.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/ndarrays.hpp>

// inversion
#include "profile_inversion.hpp"

Real ProfileInversion::LogPosteriorProbability(Radiation *prad, Hydro *phydro,
                                               Real const *par, Real *val,
                                               int k) const {
  Application::Logger app("inversion");
  int is = pmy_block_->is, ie = pmy_block_->ie;
  int ks = pmy_block_->ks;

  int ndim = GetDims();
  int nvalue = GetValues();

  // logging
  app->Log("Calculate LogPosteriorProbability");
  app->Log("I am walker = " + std::to_string(k - ks));

  Real **XpSample;
  NewCArray(XpSample, 1 + NVAPOR, plevel_.size());
  std::fill(*XpSample, *XpSample + (1 + NVAPOR) * plevel_.size(), 0.);

  // sample temperature, sample composition #1, sample composition #2, ...
  // app->Log("parameters = ", par, ndim);

  int ip = 0;
  for (auto m : idx_) {
    XpSample[m][0] = 0.;
    for (int i = 1; i <= plevel_.size() - 2; ++i)
      XpSample[m][i] = par[ip * (plevel_.size() - 2) + i - 1];
    XpSample[m][plevel_.size() - 1] = 0.;
    ip++;
  }

  // update atmosphere based on XpSample
  UpdateProfiles(phydro, XpSample, k, jl_, ju_);

  // calculate radiation for updated profiles located at j = jl_ ... ju_
  for (int j = jl_; j <= ju_; ++j) {
    app->Log("run RT for model = " + std::to_string(j));
    prad->CalRadiance(pmy_block_, 0., k, j, is, ie + 1);
  }

  // prior probability
  Real lnprior = LogPriorProbability(XpSample);

  // posterior probability
  Eigen::VectorXd misfit(nvalue);

  // calculate model result for profile at j = ju_
  CalculateFitTarget(prad, val, nvalue, k, ju_);

  Real lnpost = 0.;
  if (target_.size() > 0) {
    for (int m = 0; m < nvalue; ++m) misfit(m) = val[m] - target_(m);
    lnpost = -0.5 * misfit.transpose() * icov_ * misfit;
  }

  // posterior probability
  // Real lnprior = 0., lnpost = 0.;
  app->Log("log posterir probability = " + std::to_string(lnpost));

  FreeCArray(XpSample);

  return lnprior + lnpost;
}
