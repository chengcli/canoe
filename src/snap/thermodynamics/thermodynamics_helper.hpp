#ifndef SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HELPER_HPP_
#define SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HELPER_HPP_

// C/C++
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

// thermodynamics
#include "thermodynamics.hpp"

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

void update_gamma(Real *gamma, Real const q[]);

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

//! Potential temperature
template <typename T>
Real PotentialTemp(T w, Real p0, Thermodynamics const *pthermo) {
  Real chi = pthermo->GetChi(w);
  Real temp = pthermo->GetTemp(w);
  return temp * pow(p0 / w[IPR], chi);
}

//! Moist static energy
template <typename T>
Real MoistStaticEnergy(T w, Real gz, Thermodynamics const *pthermo) {
  Real temp = pthermo->GetTemp(w);
  Real IE = w[IDN] * pthermo->getSpecificCp(w) * temp;
  Real rho = w[IDN];
  /*if (ppart != nullptr) {
    for (int n = 0; n < NVAPOR; ++n) {
      for (int t = 0; t < ppart->u.GetDim4(); ++t) {
        rho += ppart->u(t,k,j,i);
        IE -= ppart->u(t,k,j,i)*pthermo->GetLatent(1+NVAPOR+n);
        IE += ppart->u(t,k,j,i)*ppart->GetCv(t)*temp;
      }
      ppart = ppart->next;
    }
  }*/
  return IE / rho + gz;
}

//! Relative humidity
template <typename T>
Real RelativeHumidity(T w, int iv, Thermodynamics const *pthermo) {
  Real dw[1 + NVAPOR];
  pthermo->SaturationSurplus(dw, w, VariableType::prim);
  return w[iv] / (w[iv] - dw[iv]);
}

//! Equivalent potential temperature
/*Real MoistEntropy(AthenaArray<Real> const& w, Thermodynamics *pthermo,
Particles *ppart, int k, int j, int i) { #if (NVAPOR > 0) Real gamma =
pthermo->pmy_block->peos->GetGamma(); Real tem[1] = {pthermo->GetTemp(prim)};
  update_gamma(gamma, tem);
  Real cpd = Rd_*gamma/(gamma - 1.);
  Real temp = tem[0]
  Real pres = prim[IPR];

  Real qd = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qd -= prim[n];

  Real qc[1+NVAPOR];
  std::fill(qc, qc + 1 + NVAPOR, 0.);
  for (int n = 1 + NVAPOR; n <= NVAPOR; ++n)
    qc[1+(n-1)%NVAPOR] += prim[n] + 1.0E-10;  // prevent devide by 0

  Real lv = 0.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    lv += GetLatent(n,temp)*prim[ng]*ratio;
  }

  Real st = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    st += (prim[n] + prim[ng]*ratio)*(GetCpRatio(n) - 1.);
  }
  Real lv_ov_cpt = lv/(cpd*st*temp);

  Real chi = Rd_/cpd*qd/st;

  Real xv = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    xv += prim[n]/qd/mu_ratios_[n];

  Real pd = pres/xv;

  Real rh = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    Real eta = prim[n]/qd/mu_ratios_[n];
    Real pv = pres*eta/xv;
    int nc = n + NVAPOR;
    Real esat;
    if (n == AMMONIA_VAPOR_ID)
      esat = sat_vapor_p_NH3_BriggsS(temp);
    else if (n == WATER_VAPOR_ID)
      esat = sat_vapor_p_H2O_BriggsS(temp);
    else
      esat = SatVaporPresIdeal(temp/t3_[nc], p3_[nc], beta_[nc], delta_[nc]);
    rh *= pow(pv/esat, -eta*Rd_/(cpd*st));
  }

  return temp*pow(p0/pd, chi)*exp(lv_ov_cpt)*rh;
#else
  return PotentialTemp(prim, p0, pthermo);
#endif
}*/

Real saha_ionization_electron_density(Real T, Real num, Real ion_ev);

#endif  //  SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HELPER_HPP_
