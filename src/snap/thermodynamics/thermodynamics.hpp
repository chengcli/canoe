#pragma once

// C/C++
#include <array>
#include <cfloat>
#include <iosfwd>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>

// snap
#include <snap/stride_iterator.hpp>

class MeshBlock;
class ParameterInput;

namespace Cantera {
class ThermoPhase;
class Condensation;
class Kinetics;
}  // namespace Cantera

// Thermodynamic variables are ordered in the array as the following
// 0: dry air (non-condensible gas)
// 1..NVAPOR: moist air (condensible gas)
// NVAPOR+1..NVAPOR+NCLOUD: clouds

class Thermodynamics {
  friend std::shared_ptr<Cantera::Condensation> get_kinetics_object(
      Thermodynamics const *pthermo);

 protected:
  //! Constructor for class sets up the initial conditions
  //! Protected ctor access thru static member function Instance
  Thermodynamics() {}
  static Thermodynamics *fromYAMLInput(std::string const &fname);

 public:
  enum { Size = 1 + NVAPOR + NCLOUD };

  //! thermodynamics input key in the input file [thermodynamics_config]
  static const std::string input_key;

  // member functions
  ~Thermodynamics();

  //! Return a pointer to the one and only instance of Thermodynamics
  static Thermodynamics const *GetInstance();

  static Thermodynamics const *InitFromYAMLInput(std::string const &fname);

  static Thermodynamics const *InitFromAthenaInput(ParameterInput *pin);

  //! Destroy the one and only instance of Thermodynamics
  static void Destroy();

  size_t SpeciesIndex(std::string const &name) const;

  void UpdateThermoProperties();

  //! Ideal gas constant of dry air in [J/(kg K)]
  //! \return $R_d=\hat{R}/\mu_d$
  Real GetRd() const { return Rd_; }

  //! reference adiabatic index of dry air [1]
  Real GetGammad() const { return gammad_; }

  //! mean molecular weight in [kg/mol]
  template <typename T>
  Real GetMu(T w) const;

  //! mean heat capacity in [J/(kg.K)]
  template <typename T>
  Real GetCv(T w) const;

  //! Construct an 1d atmosphere
  //! \param[in,out] qfrac mole fraction representation of air parcel
  //! \param[in] dzORdlnp vertical grid spacing
  //! \param[in] method choose from [reversible, pseudo, dry, isothermal]
  //! \param[in] grav gravitational acceleration
  //! \param[in] userp user parameter to adjust the temperature gradient
  void Extrapolate_inplace(Real dzORdlnp, std::string method, Real grav = 0.,
                           Real userp = 0.) const;

  //! Thermodnamic equilibrium at current TP
  //! \param[in,out] qfrac mole fraction representation of air parcel
  void EquilibrateTP(Real temp, Real pres) const;

  //! Thermodnamic equilibrium at current UV
  void EquilibrateUV() const;

  template <typename T>
  void SetMassFractions(StrideIterator<T *> w) const;

  template <typename T>
  void SetPrimitive(StrideIterator<T *> w) const;
  template <typename T>
  void GetPrimitive(StrideIterator<T *> w) const;

  template <typename T>
  void SetConserved(StrideIterator<T *> u, StrideIterator<T *> m) const;
  template <typename T>
  void GetConserved(StrideIterator<T *> u, StrideIterator<T *> m) const;

 public:
  Real GetTemp() const;
  Real GetPres() const;
  Real GetDensity() const;
  Real RovRd() const;

  //! \brief Inverse of the mean molecular weight
  //!
  //! Eq.16 in Li2019
  //! $ \frac{R}{R_d} = \frac{\mu_d}{\mu}$
  //! \return $1/\mu$
  template <typename T>
  Real RovRd(T w) const;

  //! \brief Effective adiabatic index
  //!
  //! Eq.71 in Li2019
  template <typename T>
  Real GetChi(T w) const;

  //! \brief Calculate temperature from primitive variable
  //!
  //! $T = p/(\rho R) = p/(\rho \frac{R}{R_d} Rd)$
  //! \return $T$
  template <typename T>
  Real GetTemp(T w) const {
    return w[IPR] / (w[IDN] * Rd_ * RovRd(w));
  }

  //! Pressure from conserved variables
  //! \return $p$
  template <typename T>
  Real GetPres(StrideIterator<T *> u, StrideIterator<T *> m) const;

  //! \brief Effective polytropic index
  //!
  //! Eq.63 in Li2019
  //! $\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  template <typename T>
  Real GetGamma(T w) const;

  template <typename T>
  Real GetEnthalpy(T w) const;

  template <typename T>
  Real GetEntropy(T w) const;

  template <typename T>
  Real GetInternalEnergy(T w) const;

  //! \brief Calculate equivalent potential temperature from primitive variable
  //!
  //! Eq.4.5.11 in Emanuel (1994)
  //! $\theta_e = T(\frac{p}{p_d})^{Rd/(cpd + cl r_t} \exp(\frac{L_v q_v}{c_p
  //! T})$
  template <typename T>
  Real MoistEntropy(T w) const;

  //! \brief Moist enthalpy
  //!
  //! Eq.66 in Li2019
  //! $h_s = c_{pd}T +  L_vq_v + L_s\sum_i q_i$
  //! \return $h_s$
  template <typename T>
  Real MoistEnthalpy(T w) const;

  template <typename T>
  std::vector<Real> SaturationSurplus(T w);

  template <typename T>
  Real IntEngToPres(StrideIterator<T *> w, Real intEng, Real rho = 1.) const {
    // internal energy
    Real fsig = 1., feps = 1.;

    // vapors
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += w[n] / rho * (cv_ratio_[n] - 1.);
      feps += w[n] / rho * (inv_mu_ratio_[n] - 1.);
    }
    // clouds
    for (int n = 1 + NVAPOR; n < Size; ++n) {
      fsig += w[n] / rho * (cv_ratio_[n] - 1.);
      feps -= w[n] / rho;
    }

    return (gammad_ - 1.) * intEng * feps / fsig;
  }

  template <typename T>
  Real PresToIntEng(StrideIterator<T *> w, Real pres, Real rho = 1.) const {
    // internal energy
    Real fsig = 1., feps = 1.;

    // vapors
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += w[n] / rho * (cv_ratio_[n] - 1.);
      feps += w[n] / rho * (inv_mu_ratio_[n] - 1.);
    }
    // clouds
    for (int n = 1 + NVAPOR; n < Size; ++n) {
      fsig += w[n] / rho * (cv_ratio_[n] - 1.);
      feps -= w[n] / rho;
    }

    return pres * fsig / feps / (gammad_ - 1.);
  }

 public:  // air parcel deprecated functions
  Real GetCvRatioMole(int n) const { return 0.; }
  Real GetCpRatioMole(int n) const { return 0.; }
  Real GetLatentEnergyMole(int n) const { return 0.; }
  Real GetCvMassRef(int n) const { return 0.; }
  Real GetInvMu(int n) const { return 0.; }

  //! Ratio of molecular weights [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\epsilon_i^{-1} =\mu_d/\mu_i$
  Real GetInvMuRatio(int n) const { return inv_mu_ratio_[n]; }

  //! Ratio of specific heat capacity [J/(kg K)] at constant volume
  //! \param[in] n the index of the thermodynamic species
  //! \return $c_{v,i}/c_{v,d}$
  Real GetCvRatio(int n) const { return cv_ratio_[n]; }

 protected:
  void _rk4_integrate_lnp(Real dlnp, std::string method, Real adlnTdlnP) const;

  void _rk4_integrate_z(Real dlnp, std::string method, Real grav,
                        Real adlnTdlnP) const;

 protected:
  std::shared_ptr<Cantera::Condensation> kinetics_;

  //! Gas constant of dry air in [J/(kg K)]
  Real Rd_;

  //! polytropic index of dry air
  Real gammad_;

  //! inverse ratio of mean molecular weights
  std::array<Real, Size> inv_mu_ratio_;

  //! ratio of specific heat capacities [J/kg] at constant pressure
  std::array<Real, Size> cp_ratio_;

  //! ratio of specific heat capacities [J/kg] at constant volume
  std::array<Real, Size> cv_ratio_;

  //! pointer to the single Thermodynamics instance
  static Thermodynamics *mythermo_;
};

//! \brief Calculate moist adiabatic temperature gradient
//!
//! Eq.68 in Li2019
//! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
//! \return $\Gamma_m$
Real cal_dlnT_dlnP(Real const *xfrac, Real gammad, Real const *cp_ratio_mole,
                   Real const *latent);

#include "thermodynamics_athena.hpp"
