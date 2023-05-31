// C/C++
#include <athena/athena_arrays.hpp>
#include <athena/hydro/hydro.hpp>
#include <cstdlib>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"
#include "thermodynamics_helper.hpp"

void Thermodynamics::ConstructAtmosphere(Real **w, Real Ts, Real Ps, Real grav,
                                         Real dzORdlnp, int len, Adiabat method,
                                         Real userp) const {
  Real gamma = pmy_block->peos->GetGamma();

  // hydro + vapor + condensate
  Real q1[NHYDRO + 2 * NVAPOR];
  w[0][IDN] = w[0][IPR] = 1.;
  PrimitiveToChemical(q1, w[0]);
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) q1[n] = 0.;
  // change molar concentration to mixing ratio
  Real mols = q1[IPR] / (q1[IDN] * Rgas);
  for (int n = 1; n <= NVAPOR; ++n) q1[n] /= mols;
  // reset TP
  q1[IDN] = Ts;
  q1[IPR] = Ps;

  Real rcp[1 + 3 * NVAPOR];  // molar cp ratio
  for (int n = 0; n <= 3 * NVAPOR; ++n) rcp[n] = cp_ratios_[n] * mu_ratios_[n];

  int isat[1 + NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n) isat[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  Real xg = 1.;  // total moles in gas phase
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    int nc = q1[IDN] > t3_[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
    int ic = NHYDRO - NVAPOR + nc - 1;
    Real rate = VaporCloudEquilibrium(q1, iv, ic, t3_[iv], p3_[iv], 0.,
                                      beta_[nc], delta_[nc]);
    q1[iv] -= rate;
    q1[ic] += rate;
    if (rate > 0.) isat[iv] = 1;
    xg -= rate;
  }

  // change molar mixing ratio to molar concentration
  Real qv = 1.;
  for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) qv -= q1[n];
  mols = q1[IPR] / (q1[IDN] * Rgas) / qv;
  for (int n = 1; n <= NVAPOR; ++n) q1[n] *= mols;
  // set vapor
  ChemicalToPrimitive(w[0], q1);
  // change back
  for (int n = 1; n <= NVAPOR; ++n) q1[n] /= mols;
  // set clouds, mass density
  Real mu_d = Rgas / Rd_;
  for (int n = 0; n < 2 * NVAPOR; ++n)
    w[0][NHYDRO + n] =
        q1[NHYDRO + n] * mols * mu_ratios_[1 + NVAPOR + n] * mu_d;

  for (int i = 1; i < len; ++i) {
    // RK4 integration
    Real rdlnTdlnP, adTdz;
#ifdef HYDROSTATIC
    rdlnTdlnP = userp;
    rk4_integrate_lnp_adaptive(q1, isat, rcp, beta_, delta_, t3_, p3_, gamma,
                               dzORdlnp, ftol_, static_cast<int>(method),
                               rdlnTdlnP);
#else
    adTdz = userp;
    rk4_integrate_z_adaptive(q1, isat, rcp, mu_ratios_, beta_, delta_, t3_, p3_,
                             gamma, grav / Rd_, dzORdlnp, ftol_,
                             static_cast<int>(method), adTdz);
#endif
    // reset mols
    qv = 1.;
    for (int n = NHYDRO; n < NHYDRO + 2 * NVAPOR; ++n) qv -= q1[n];
    mols = q1[IPR] / (q1[IDN] * Rgas) / qv;
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
}
