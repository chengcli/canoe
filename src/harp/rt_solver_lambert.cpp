// C/C++
#include <cmath>
#include <iostream>

// athena
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/special.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// harp
#include "radiation.hpp"
#include "rt_solvers.hpp"

void RadiationBand::RTSolverLambert::CalBandFlux(MeshBlock const *pmb, int k,
                                                 int j, int il, int iu) {
  throw NotImplementedError("RTSolverLambert::CalBandFlux");
}

void RadiationBand::RTSolverLambert::CalBandRadiance(MeshBlock const *pmb,
                                                     int k, int j) {
  RadiationBand *pband = pmy_band_;

  auto &rayOutput = pband->rayOutput_;
  auto &toa = pband->toa_;
  auto &tau = pband->tau_;
  auto &temf = pband->temf_;
  auto &spec = pband->pgrid_->spec;

  auto &btoa = pband->btoa;
  //! \note $T ~ Ts*(\tau/\tau_s)^\alpha$ at lower boundary
  Real alpha = pband->HasPar("alpha") ? pband->GetPar<Real>("alpha") : 0.;

  Real il = pmb->is;
  Real iu = pmb->ie + 1;

  std::vector<Real> taut(iu + 1);

  // integrate from top to bottom
  for (int m = 0; m < pband->GetNumOutgoingRays(); ++m) {
    btoa(m, k, j) = 0.;
    for (int n = 0; n < pband->GetNumSpecGrids(); ++n) {
      taut[iu] = 0.;
      toa(n, m) = 0.;
      for (int i = iu - 1; i >= il; --i) {
        taut[i] = taut[i + 1] + tau(n, i) / rayOutput[m].mu;
        toa(n, m) +=
            0.5 * (temf[i + 1] * exp(-taut[i + 1]) + temf[i] * exp(-taut[i])) *
            tau(n, i) / rayOutput[m].mu;
      }
      toa(n, m) += temf[il] * exp(-taut[il]);
      //! \note correction for small optical opacity
      if ((alpha > 0) && (taut[il] < 1000.)) {
        toa(n, m) += temf[il] * alpha * gammq(alpha, taut[il]) *
                     pow(taut[il], -alpha) * tgamma(alpha);
      }
      btoa(m, k, j) += spec[n].wght * toa(n, m);
    }
  }
}
