// C/C++
#include <iostream>

// canoe
#include <configure.hpp>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/meshblock_impl.hpp>

// debugger
#include <debugger/debugger.hpp>

// inversion
#include "gaussian_process.hpp"
#include "profile_inversion.hpp"

extern std::unique_ptr<Debugger> pdebug;

Real ProfileInversion::LogPriorProbability(Real **XpSample) const {
  int nsample = plevel_.size();
  std::vector<Real> zlev(nsample);
  std::vector<Real> zstd(nsample);

  Real P0 = pmy_block_->pimpl->GetReferencePressure();
  Real H0 = pmy_block_->pimpl->GetPressureScaleHeight();

  for (int i = 0; i < nsample; ++i) zlev[i] = -H0 * log(plevel_[i] / P0);

  Real lnprior = 0.;
  for (auto m : idx_) {
    for (int i = 0; i < nsample; ++i)
      zstd[i] = Xstd_[m] * pow(exp(zlev[i] / H0), chi_);
    lnprior += gp_lnprior(SquaredExponential, XpSample[m], zlev.data(),
                          zstd.data(), nsample, Xlen_[m]);
  }

  pdebug->Message("Log prior probability", lnprior);

  return lnprior;
}
