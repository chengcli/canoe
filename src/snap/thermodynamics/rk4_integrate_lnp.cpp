// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"
#include "thermodynamics_helper.hpp"

void rk4_integrate_lnp(Real q[], int isat[], Real const rcp[],
                       Real const beta[], Real const delta[], Real const t3[],
                       Real const p3[], Real gamma, Real dlnp, int method,
                       Real rdlnTdlnP) {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = q[IDN];
  Real pres = q[IPR];
  Real chi[4];

  Real beta0[1 + 3 * NVAPOR], delta0[1 + 3 * NVAPOR];
  std::fill(beta0, beta0 + 1 + 3 * NVAPOR, 0.);
  std::fill(delta0, delta0 + 1 + 3 * NVAPOR, 0.);

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor and cloud
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      q[iv] += q[NHYDRO + iv - 1] + q[NHYDRO + NVAPOR + iv - 1];
      q[NHYDRO + iv - 1] = 0.;
      q[NHYDRO + NVAPOR + iv - 1] = 0.;
      isat[iv] = 0;
    }

    for (int iv = 1; iv <= NVAPOR; ++iv) {
      int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
      int ic = NHYDRO - NVAPOR + nc - 1;
      Real rate = VaporCloudEquilibrium(q, iv, ic, t3[iv], p3[iv], 0., beta[nc],
                                        delta[nc]);
      if (rate > 0.) isat[iv] = 1;
      q[iv] -= rate;
      if (method == 0) q[ic] += rate;
    }

    // calculate gamma
    update_gamma(&gamma, q);

    // calculate tendency
    if (method == 0 || method == 1)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta, delta, t3, gamma);
    else if (method == 2)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta0, delta0, t3, gamma);
    else  // isothermal
      chi[rk] = 0.;
    chi[rk] += rdlnTdlnP - 1.;

    // integrate over dlnp
    if (rk < 3)
      q[IDN] = temp * exp(chi[rk] * dlnp * step[rk]);
    else
      q[IDN] = temp * exp(1. / 6. *
                          (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp);
    q[IPR] = pres * exp(dlnp);
  }

  // recondensation
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
    int ic = NHYDRO - NVAPOR + nc - 1;
    Real rate = VaporCloudEquilibrium(q, iv, ic, t3[iv], p3[iv], 0., beta[nc],
                                      delta[nc]);
    if (rate > 0.) isat[iv] = 1;
    q[iv] -= rate;
    if (method == 0) q[ic] += rate;
  }
}

void rk4_integrate_lnp_adaptive(Real q[], int isat[], Real const rcp[],
                                Real const beta[], Real const delta[],
                                Real const t3[], Real const p3[], Real gamma,
                                Real dlnp, Real ftol, int method,
                                Real rdlnTdlnP) {
  Real q1[NHYDRO + 2 * NVAPOR], q2[NHYDRO + 2 * NVAPOR];
  int isat1[1 + NVAPOR], isat2[1 + NVAPOR];

  for (int n = 0; n < NHYDRO + 2 * NVAPOR; ++n) {
    q1[n] = q[n];
    q2[n] = q[n];
  }
  for (int n = 0; n <= NVAPOR; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_lnp(q1, isat1, rcp, beta, delta, t3, p3, gamma, dlnp, method,
                    rdlnTdlnP);

  // refined step
  rk4_integrate_lnp(q2, isat2, rcp, beta, delta, t3, p3, gamma, dlnp / 2.,
                    method, rdlnTdlnP);
  rk4_integrate_lnp(q2, isat2, rcp, beta, delta, t3, p3, gamma, dlnp / 2.,
                    method, rdlnTdlnP);

  if (fabs(q2[IDN] - q1[IDN]) > ftol) {
    rk4_integrate_lnp_adaptive(q, isat, rcp, beta, delta, t3, p3, gamma,
                               dlnp / 2., ftol / 2., method, rdlnTdlnP);
    rk4_integrate_lnp_adaptive(q, isat, rcp, beta, delta, t3, p3, gamma,
                               dlnp / 2., ftol / 2., method, rdlnTdlnP);
  } else {
    for (int n = 0; n < NHYDRO + 2 * NVAPOR; ++n) q[n] = q2[n];
    for (int n = 0; n <= NVAPOR; ++n) isat[n] = isat2[n];
  }
}
