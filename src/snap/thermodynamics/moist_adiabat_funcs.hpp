#ifndef MOIST_ADIABAT_FUNCS_HPP
#define MOIST_ADIABAT_FUNCS_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include <athena/athena.hpp>

/*! Ideal saturation vapor pressure
 *
 * $p^* = p^r\exp[\beta(1-1/t)-\delta\lnt]$
 * $p^r$ is the reference pressure, usually choosen to be the triple point
 * pressure \n $t=T/T^r$ is the dimensionless temperature. $T^r$ is the
 * reference temperature. \n Similar to $p^r$, $T^r$ is usually choosen to be
 * the triple point temperature
 * @return $p^*$ [pa]
 */
inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real delta) {
  return p3 * exp((1. - 1. / t) * beta - delta * log(t));
}

// Real r_ov_rd(Real const w[], Real const eps[]) {
//   Real feps = 1.;
//   for (int n = 1; n <= NVAPOR; ++n)
//     feps += w[n]*(1./eps[n] - 1.);
//   return feps;
// }
// Real q_gas(Real const w[]);
// Real q_sig(Real const w[], Real const sig[]);
// Real qhat_eps(Real const q[], Real const eps[]);
// Real qhat_rcp(Real const q[], Real const rcp[]);

/*! Change mass mixing ratio to molar mixing ratio
inline void mass_to_molar(Real q[], Real const w[], Real const eps[], Real Rd) {
  // set molar mixing ratio
  Real sum = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    q[n] = w[n]/eps[n];
    sum += w[n]*(1./eps[n] - 1.);
  }

  for (int n = 1; n <= NVAPOR; ++n)
    q[n] /= sum;

  // set pressure and temperature
  q[IPR] = w[IPR];
  q[IDN] = w[IPR]/(w[IDN]*Rd*sum);
}

 Change molar mixing ratio to mass mixing ratio
inline void molar_to_mass(Real w[], Real const q[], Real const eps[], Real Rd) {
  // set mass mixing ratio
  Real sum = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] = q[n]*eps[n];
    sum += q[n]*(eps[n] - 1.);
  }
  for (int n = 1; n <= NVAPOR; ++n)
    w[n] /= sum;

  // set pressure and density
  w[IPR] = q[IPR];
  w[IDN] = sum*q[IPR]/(q[IDN]*Rd);
}*/

/*! Moist adiabatic lapse rate
 *
 * Calculate moist adiabatic lapse rate according to @cite li2018
 * @param q an array of [T,q1,q2,...,qn,P,q11,q12,q21,q22,...,qn1,qn2]
 * @param isat indicator array of saturated or not
 * @param rcp molar cp ratios with respect to dry air
 */
Real dlnTdlnP(Real const q[], int const isat[], Real const rcp[],
              Real const beta[], Real const delta[], Real const t3[],
              Real gamma);

void update_gamma(Real &gamma, Real const q[]);

/*! Calculate the equilibrium between vapor and cloud
 *
 * @param q molar mixing ratio
 * #param iv index of vapor
 * @param ic index of cloud
 * @param t3 triple point temperature
 * @param p3 triple point pressure
 * @param alpha = L/cv evaluated at current temperature
 * @return molar change of vapor to cloud
 */
Real VaporCloudEquilibrium(Real const q[], int iv, int ic, Real t3, Real p3,
                           Real alpha, Real beta, Real delta,
                           bool no_cloud = false);

void rk4_integrate_z(Real q[], int isat[], Real rcp[], Real const eps[],
                     Real const beta[], Real const delta[], Real const t3[],
                     Real const p3[], Real gamma, Real g_ov_Rd, Real dz,
                     int method, Real adTdz = 0.);

void rk4_integrate_z_adaptive(Real q[], int isat[], Real rcp[],
                              Real const eps[], Real const beta[],
                              Real const delta[], Real const t3[],
                              Real const p3[], Real gamma, Real g_ov_Rd,
                              Real dz, Real ftol, int method, Real adTdz = 0.);

void rk4_integrate_lnp(Real q[], int isat[], Real const rcp[],
                       Real const beta[], Real const delta[], Real const t3[],
                       Real const p3[], Real gamma, Real dlnp, int method,
                       Real rdlnTdlnP = 1.);

void rk4_integrate_lnp_adaptive(Real q[], int isat[], Real const rcp[],
                                Real const beta[], Real const delta[],
                                Real const t3[], Real const p3[], Real gamma,
                                Real dlnp, Real ftol, int method,
                                Real rdlnTdlnP = 1.);

bool read_idealized_parameters(std::string fname, Real &pmin, Real &pmax,
                               AthenaArray<Real> &dTdz, Real *&pmin_q,
                               Real *&pmax_q, AthenaArray<Real> *&dqdz,
                               std::vector<int> &mindex);

#endif
