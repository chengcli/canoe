#include <algorithm>
#include <iostream>

#include "moist_adiabat_funcs.hpp"

void rk4_integrate_z(Real q[], int isat[], Real rcp[], Real const eps[],
                     Real const beta[], Real const delta[], Real const t3[],
                     Real const p3[], Real gamma, Real g_ov_Rd, Real dz,
                     int method, Real adTdz) {
  Real step[] = {0.5, 0.5, 1.};
  Real temp = q[IDN];
  Real pres = q[IPR];
  Real dTdz[4], chi[4];

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

    Real q_gas = 1., q_eps = 1.;
    for (int iv = 1; iv <= NVAPOR; ++iv) {
      int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
      int ic = NHYDRO - NVAPOR + nc - 1;
      Real rate = VaporCloudEquilibrium(q, iv, ic, t3[iv], p3[iv], 0., beta[nc],
                                        delta[nc]);
      // Real rate = 0.;
      if (rate > 0.) isat[iv] = 1;
      q[iv] -= rate;
      if (method == 0) q[ic] += rate;
      q_gas -= q[ic];
      q_eps += (q[iv] + q[ic]) * (eps[iv] - 1.);
    }
    Real R_ov_Rd = q_gas / q_eps;

    // calculate gamma
    update_gamma(gamma, q);

    if (method == 0 || method == 1)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta, delta, t3, gamma);
    else if (method == 2)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta0, delta0, t3, gamma);
    else  // isothermal
      chi[rk] = 0.;
    dTdz[rk] = -chi[rk] * g_ov_Rd / R_ov_Rd + adTdz;
    chi[rk] = -R_ov_Rd / g_ov_Rd * dTdz[rk];

    // integrate over dz
    Real chi_avg;
    if (rk < 3) {
      q[IDN] = temp + dTdz[rk] * dz * step[rk];
      chi_avg = chi[rk];
    } else {
      q[IDN] = temp +
               1. / 6. * (dTdz[0] + 2. * dTdz[1] + 2. * dTdz[2] + dTdz[3]) * dz;
      chi_avg = 1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]);
    }
    if (!(q[IDN] > 0.)) q[IDN] = temp;
    if (fabs(q[IDN] - temp) > 0.1)  // isothermal limit
      q[IPR] = pres * pow(q[IDN] / temp, 1. / chi_avg);
    else
      q[IPR] = pres * exp(-2. * g_ov_Rd * dz / (R_ov_Rd * (q[IDN] + temp)));
  }

  // recondensation
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
    int ic = NHYDRO - NVAPOR + nc - 1;
    Real rate = VaporCloudEquilibrium(q, iv, ic, t3[iv], p3[iv], 0., beta[nc],
                                      delta[nc]);
    // Real rate = 0.;
    if (rate > 0.) isat[iv] = 1;
    q[iv] -= rate;
    if (method == 0) q[ic] += rate;
  }
}

void rk4_integrate_z_adaptive(Real q[], int isat[], Real rcp[],
                              Real const eps[], Real const beta[],
                              Real const delta[], Real const t3[],
                              Real const p3[], Real gamma, Real g_ov_Rd,
                              Real dz, Real ftol, int method, Real adTdz) {
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
  rk4_integrate_z(q1, isat1, rcp, eps, beta, delta, t3, p3, gamma, g_ov_Rd, dz,
                  method, adTdz);

  // refined step
  rk4_integrate_z(q2, isat2, rcp, eps, beta, delta, t3, p3, gamma, g_ov_Rd,
                  dz / 2., method, adTdz);
  rk4_integrate_z(q2, isat2, rcp, eps, beta, delta, t3, p3, gamma, g_ov_Rd,
                  dz / 2., method, adTdz);

  // abort if dz is less than 0.001 m
  // std::cout << dz << " ";
  // std::cout << fabs(q2[IDN] - q1[IDN]) << " ";
  // std::cout << ftol << std::endl;
  // assert(dz > 0.001);

  if (fabs(q2[IDN] - q1[IDN]) > ftol) {
    rk4_integrate_z_adaptive(q, isat, rcp, eps, beta, delta, t3, p3, gamma,
                             g_ov_Rd, dz / 2., ftol, method, adTdz);
    rk4_integrate_z_adaptive(q, isat, rcp, eps, beta, delta, t3, p3, gamma,
                             g_ov_Rd, dz / 2., ftol, method, adTdz);
  } else {
    for (int n = 0; n < NHYDRO + 2 * NVAPOR; ++n) q[n] = q2[n];
    for (int n = 0; n <= NVAPOR; ++n) isat[n] = isat2[n];
  }
}
