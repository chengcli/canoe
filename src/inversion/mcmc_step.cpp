// C/C++
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

// canoe
#include <configure.hpp>

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

void Inversion::MCMCMove(Radiation *prad, Hydro *phydro) {
  Application::Logger app("inversion");
  app->Log("MCMC move");

  // initialize model
  if (recs_.cur == 0) {
    MCMCInit(prad, phydro);
  }

  std::stringstream msg;
  if (!mcmc_initialized_) {
    app->Error("mcmc chain uninitialized");
  }

  int cur = recs_.cur;
  int ndim = recs_.ndim;
  int nwalker = recs_.nwalker;

  for (int k = 0; k < nwalker; ++k) {
    zz_[k] =
        mcmc_stretch_move(par_, recs_.par[cur - 1], k, nwalker, ndim, &opts_);
    // mcmc_walk_move(par_all + k*np, par, np, k, nwalker, zz + k, opts_);

    recs_.lnp[cur][k] = LogPosteriorProbability(
        prad, phydro, par_, recs_.val[cur][k], pmy_block_->ks + k);
  }
}

void Inversion::MCMCSave(Hydro *phydro) {
  int cur = recs_.cur;
  int ndim = recs_.ndim;
  int nwalker = recs_.nwalker;
  int nval = recs_.nvalue;

  int is = pmy_block_->is, ie = pmy_block_->ie;
  int ks = pmy_block_->ks;

  unsigned int seed = time(NULL);
  for (int k = 0; k < nwalker; ++k) {
    double lnp0 = recs_.lnp[cur - 1][k], lnp1 = recs_.lnp[cur][k],
           pdiff = (ndim - 1.) * log(zz_[k]) + lnp1 - lnp0;

    if (pdiff > log(1. * rand_r(&seed) / RAND_MAX)) {  // accept this position
      for (int d = 0; d < ndim; ++d) recs_.par[cur][k][d] = par_[d];

      if (lnp1 > recs_.opt_lnp) {  // new best position
        recs_.opt_lnp = lnp1;
        for (int d = 0; d < ndim; ++d) recs_.opt_par[d] = recs_.par[cur][k][d];
        for (int d = 0; d < nval; ++d) recs_.opt_val[d] = recs_.val[cur][k][d];
      }

      recs_.accept++;
      recs_.newstate[cur][k] = 1;

      // save hydro to w1
      for (int n = 0; n < NHYDRO; ++n)
        for (int j = jl_; j <= ju_; ++j)
          for (int i = is; i <= ie; ++i) {
            phydro->w1(n, ks + k, j, i) = phydro->w(n, ks + k, j, i);
          }
    } else {  // do not accept this position
      for (int d = 0; d < ndim; ++d)
        recs_.par[cur][k][d] = recs_.par[cur - 1][k][d];
      for (int d = 0; d < nval; ++d)
        recs_.val[cur][k][d] = recs_.val[cur - 1][k][d];
      recs_.lnp[cur][k] = recs_.lnp[cur - 1][k];
      recs_.newstate[cur][k] = 0;

      // load hydro from w1
      for (int n = 0; n < NHYDRO; ++n)
        for (int j = jl_; j <= ju_; ++j)
          for (int i = is; i <= ie; ++i) {
            phydro->w(n, ks + k, j, i) = phydro->w1(n, ks + k, j, i);
          }
    }
  }

  recs_.cur++;
  if (recs_.cur % opts_.print == 0) mcmc_report(&opts_, &recs_, "a");
}
