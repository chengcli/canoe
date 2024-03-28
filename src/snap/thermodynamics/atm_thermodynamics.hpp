#ifndef SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_
#define SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// snap
#include "thermodynamics.hpp"

//! adiaibatic index of dry air [1]
Real get_gammad(AirParcel const &qfrac);

//! Adiabatic index
//! \return $\chi = \frac{R}{c_p}$
Real get_chi(AirParcel const &qfrac);

//! Inverse of the mean molecular weight (with cloud)
//! \param[in] qfrac mole fraction representation of air parcel
Real get_rovrd(AirParcel const &qfrac);

//! \brief Calculate moist adiabatic temperature gradient
//!
//! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
//! \return $\Gamma_m$
Real cal_dlnT_dlnP(Thermodynamics const *pthermo, AirParcel const &qfrac,
                   Real latent[])

    // thermodynamic functions
    //! \brief Relative humidity
    //!
    //! $H_i = \frac{e_i}{e_i^s}$
    //! \return $H$
    inline Real get_relative_humidity(Thermodynamics const *pthermo,
                                      AirParcel const &qfrac, int n) {
  qfrac.ToMoleFraction();
  auto rates = TryEquilibriumTP_VaporCloud(qfrac, n, 0., true);
  return qfrac.w[n] / (qfrac.w[n] + rates[0]);
}

//! Specific heat capacity [J/(kg K)] at constant volume
//! $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
//! $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
//! \return $c_{v,i}$
inline Real get_cv_mass(Thermodynamics const *pthermo, AirParcel const &qfrac,
                        int n) const {
  qfrac.ToMoleFraction();
  Real cvd = pthermo->GetRd() / (get_gammad(qfrac) - 1.);
  return pthermo->GetCvRatioMass(n) * cvd;
}

//! Specific heat capacity [J/(mol K)] of the air parcel at constant volume
//! \return $\hat{c}_v$
inline Real get_cv_mole(Thermodynamics const *pthermo, AirParcel const &qfrac,
                        int n) {
  qfrac.ToMoleFraction();
  return get_cv_mass(pthermo, qfrac, n) * pthermo->GetMu(n);
}

//! Specific heat capacity [J/(kg K)] at constant pressure
//! $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
//! $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
//! \return $c_p$
inline Real get_cp_mass(Thermodynamics const *pthermo, AirParcel const &qfrac,
                        int n) {
  Real gammad = get_gammad(qfrac);
  Real cpd = pthermo->GetRd() * gammad / (gammad - 1.);
  return pthermo->GetCpRatioMass(n) * cpd;
}

//! Specific heat capacity [J/(mol K)] of the air parcel at constant pressure
//! \return $\hat{c}_v$
inline Real get_cp_mole(Thermodynamics const *pthermo, AirParcel const &qfrac,
                        int n) {
  return get_cp_mass(qfrac, n) * pthermo->GetMu(n);
}

#endif  // SRC_SNAP_THERMODYNAMICS_ATM_THERMODYNAMICS_HPP_
