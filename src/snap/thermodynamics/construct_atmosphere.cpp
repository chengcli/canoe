// C/C++
#include <cstdlib>
#include <iostream>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>

// thermodynamics
#include "thermodynamics.hpp"
#include "thermodynamics_helper.hpp"

void Thermodynamics::ConstructAtmosphere(Variable *qfrac, Real dzORdlnp,
                                         Adiabat method, Real userp) const {
  Real gamma = pmy_block_->peos->GetGamma();
  Real grav = -pmy_block_->phydro->hsrc.GetG1();
  auto pimpl = pmy_block_->pimpl;

  int isat[1 + NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n) isat[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  Real xg = 1.;  // total moles in gas phase
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    int nc = qfrac->w[IDN] > t3_[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
    int ic = NHYDRO - NVAPOR + nc - 1;
    Real rate = VaporCloudEquilibrium(qfrac, iv, ic, t3_[iv], p3_[iv], 0.,
                                      beta_[nc], delta_[nc]);
    q1->w[iv] -= rate;
    q1->w[ic] += rate;
    if (rate > 0.) isat[iv] = 1;
    xg -= rate;
  }

  // RK4 integration
  Real adlnTdlnP, adTdz;
#ifdef HYDROSTATIC
  adlnTdlnP = userp;
  rk4_integrate_lnp_adaptive(q1, isat, rcp, beta_, delta_, t3_, p3_, gamma,
                             dzORdlnp, ftol_, static_cast<int>(method),
                             adlnTdlnP);
#else
  adTdz = userp;
  rk4_integrate_z_adaptive(q1, isat, rcp, mu_ratios_, beta_, delta_, t3_, p3_,
                           gamma, grav / Rd_, dzORdlnp, ftol_,
                           static_cast<int>(method), adTdz);
#endif
  // reset mols
  qv = 1.;
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) qv -= q1[n];
  mols = q1[IPR] / (q1[IDN] * Constants::Rgas) / qv;
  // change molar mixing ratio to molar concentration
  for (int n = 1; n <= NVAPOR; ++n) q1[n] *= mols;
  // set vapor
  ChemicalToPrimitive(w[i], q1);
  // change back
  for (int n = 1; n <= NVAPOR; ++n) q1[n] /= mols;
  // set clouds, mass density
  for (int n = 0; n < 2 * NVAPOR; ++n)
    w[i][NHYDRO + n] =
        q1[NHYDRO + n] * mols * mu_ratios_[1 + NVAPOR + n] * mu_d;
}
