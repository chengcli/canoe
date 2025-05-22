// C/C++ headers
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

// canoe
#include <configure.h>

// climath
#include <climath/core.h>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>

// inversion
#include "inversion.hpp"
#include "mcmc.hpp"

void Inversion::MCMCInit(Radiation *prad, Hydro *phydro) {
  int nwalker = recs_.nwalker;
  int ndim = recs_.ndim;
  int nval = recs_.nvalue;

  Application::Logger app("inversion");

  Real **par = init_pos_;

  recs_.lnp[0][0] = 1.;

  int is = pmy_block_->is, ie = pmy_block_->ie;
  int ks = pmy_block_->ks;

  // make sure that the start points are valid
  for (int k = 0; k < nwalker; ++k) {
    recs_.lnp[0][k] =
        LogPosteriorProbability(prad, phydro, par[k], recs_.val[0][k], ks + k);

    int niter = 0;
    while (std::isnan(recs_.lnp[0][k]) &&
           (niter++ < 10)) {  // if point (k) is invalid
      mcmc_stretch_move(par[k], par, k, nwalker, ndim, &opts_);
      recs_.lnp[0][k] = LogPosteriorProbability(prad, phydro, par[k],
                                                recs_.val[0][k], ks + k);
    }

    if (niter >= 10) {
      app->Error("Starting point iteration > 10 times");
    }

    // transfer input parameters to records
    for (int d = 0; d < ndim; ++d) recs_.par[0][k][d] = par[k][d];

    if (recs_.lnp[0][k] > recs_.opt_lnp) {
      recs_.opt_lnp = recs_.lnp[0][k];
      for (int d = 0; d < ndim; ++d) recs_.opt_par[d] = recs_.par[0][k][d];
      for (int d = 0; d < nval; ++d) recs_.opt_val[d] = recs_.val[0][k][d];
    }
    recs_.newstate[0][k] = 1;

    // save hydro to w1
    for (int n = 0; n < NHYDRO; ++n)
      for (int j = jl_; j <= ju_; ++j)
        for (int i = is; i <= ie; ++i) {
          phydro->w1(n, ks + k, j, i) = phydro->w(n, ks + k, j, i);
        }
  }

  recs_.cur++;
  recs_.accept += nwalker;

  // make initial output
  mcmc_report(&opts_, &recs_, "w");
}
