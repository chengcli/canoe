#ifndef SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_
#define SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_

// canoe
#include <configure.h>

#include <air_parcel.hpp>

// snap
#include "thermodynamics.hpp"

//! adiaibatic index of dry air [1]
Real get_gammad(AirParcel const &qfrac);

//! Adiabatic index
//! \return $\chi = \frac{R}{c_p}$
Real get_chi(AirParcel const &qfrac, Real const *cp_ratio_mole);

//! Inverse of the mean molecular weight (with cloud)
//! \param[in] qfrac mole fraction representation of air parcel
Real get_rovrd(AirParcel const &qfrac, Real const *mu_ratio);

void set_total_equivalent_vapor(
    AirParcel *qfrac, IndexSet const *cloud_index_set,
    std::map<IndexPair, ReactionInfo> const &cloud_reaction_map);

//! update T/P
void update_TP_conserving_U(AirParcel *qfrac, Real rmole, Real umole,
                            Real const *cv_ratio_mole,
                            Real const *latent_energy_mole,
                            IndexSet const *cloud_index_set);

Real get_internal_energy_mole(AirParcel const &qfrac, Real const *cv_ratio_mole,
                              Real const *latent_energy_mole,
                              IndexSet const *cloud_index_set);

void rk4_integrate_lnp(AirParcel *qfrac, Real dlnp, std::string method,
                       Real adlnTdlnP);

void rk4_integrate_z(AirParcel *qfrac, Real dlnp, std::string method, Real grav,
                     Real adlnTdlnP);

//! \brief Calculate moist adiabatic temperature gradient
//!
//! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
//! \return $\Gamma_m$
Real cal_dlnT_dlnP(AirParcel const &qfrac, Real const *cp_ratio_mole,
                   Real const *latent);

// thermodynamic functions
//! \brief Relative humidity
//!
//! $H_i = \frac{e_i}{e_i^s}$
//! \return $H$
inline Real get_relative_humidity(AirParcel const &qfrac, int n) {
  auto pthermo = Thermodynamics::GetInstance();
  auto rates = pthermo->TryEquilibriumTP_VaporCloud(qfrac, n, 0., true);
  return qfrac.w[n] / (qfrac.w[n] + rates[0]);
}

//! Specific heat capacity [J/(kg K)] at constant volume
//! $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
//! $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
//! \return $c_{v,i}$
inline Real get_cv_mass(AirParcel const &qfrac, int n, Real Rd,
                        Real const *cv_ratio_mass) {
  Real cvd = Rd / (get_gammad(qfrac) - 1.);
  return cv_ratio_mass[n] * cvd;
}

//! Specific heat capacity [J/(mol K)] of the air parcel at constant volume
//! \return $\hat{c}_v$
inline Real get_cv_mole(AirParcel const &qfrac, int n, Real Rd,
                        Real const *cv_ratio_mole) {
  Real cvd = Rd / (get_gammad(qfrac) - 1.);
  return cv_ratio_mole[n] * cvd;
}

//! Specific heat capacity [J/(kg K)] at constant pressure
//! $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
//! $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
//! \return $c_p$
inline Real get_cp_mass(AirParcel const &qfrac, int n, Real Rd,
                        Real const *cp_ratio_mass) {
  Real gammad = get_gammad(qfrac);
  Real cpd = Rd * gammad / (gammad - 1.);
  return cp_ratio_mass[n] * cpd;
}

//! Specific heat capacity [J/(mol K)] of the air parcel at constant pressure
//! \return $\hat{c}_v$
inline Real get_cp_mole(AirParcel const &qfrac, int n, Real Rd,
                        Real const *cp_ratio_mole) {
  Real gammad = get_gammad(qfrac);
  Real cpd = Rd * gammad / (gammad - 1.);
  return cp_ratio_mole[n] * cpd;
}

inline Real get_density_mole(AirParcel const &qfrac) {
  Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
  for (int n = 0; n < NCLOUD; ++n) qgas += -qfrac.c[n];
  return qfrac.w[IPR] / (Constants::Rgas * qfrac.w[IDN] * qgas);
}

#endif  // SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_
