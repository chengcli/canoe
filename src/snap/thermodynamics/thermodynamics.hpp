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

  //! Ratio of molecular weights [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\epsilon_i^{-1} =\mu_d/\mu_i$
  Real GetInvMuRatio(int n) const { return inv_mu_ratio_[n]; }

  //! Ratio of specific heat capacity [J/(kg K)] at constant volume
  //! \param[in] n the index of the thermodynamic species
  //! \return $c_{v,i}/c_{v,d}$
  Real GetCvRatio(int n) const { return cv_ratio_[n]; }

  //! Ratio of specific heat capacity [J/(kg K)] at constant pressure
  //! \param[in] n the index of the thermodynamic species
  //! \return $c_{p,i}/c_{p,d}$
  Real GetCpRatio(int n) const { return cp_ratio_[n]; }

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
  void EquilibrateTP() const;

  //! Thermodnamic equilibrium at current UV
  void EquilibrateUV() const;

  void SetTemperature(Real temp) const;
  void SetPressure(Real pres) const;
  void SetDensity(Real dens) const;
  template <typename T>
  void SetMassFractions(StrideIterator<T *> w) const;

  template <typename T>
  void SetPrimitive(StrideIterator<T *> w) const;
  template <typename T>
  void GetPrimitive(StrideIterator<T *> w) const;

  template <typename T>
  void SetConserved(StrideIterator<T *> u, StrideIterator<T *> m) const;
  template <typename T>
  void GetConserved(StrideIterator<T *> u) const;

 public:
  //! \brief Inverse of the mean molecular weight
  //!
  //! Eq.16 in Li2019
  //! $ \frac{R}{R_d} = \frac{\mu_d}{\mu}$
  //! \return $1/\mu$
  template <typename T>
  Real RovRd(T w) const;

  Real RovRd() const;

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

  Real GetTemp() const;

  //! Pressure from conserved variables
  //! \return $p$
  template <typename T>
  Real GetPres(StrideIterator<T *> u, StrideIterator<T *> m) const;

  Real GetPres() const;

  Real GetDensity() const;

  //! \brief Calculate potential temperature from primitive variable
  //!
  //! $\theta = T(\frac{p_0}{p})^{\chi}$
  //! \return $\theta$
  template <typename T>
  Real PotentialTemp(T w, Real p0) const {
    return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
  }

  //! \brief Effective polytropic index
  //!
  //! Eq.63 in Li2019
  //! $\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  template <typename T>
  Real GetGamma(T w) const;

  template <typename T>
  Real GetEntropy(T w) const;

  template <typename T>
  Real GetEnthalpy(T w) const;

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

 public:  // air parcel deprecated functions
  Real GetCvRatioMole(int n) const { return 0.; }
  Real GetCpRatioMole(int n) const { return 0.; }
  Real GetLatentEnergyMole(int n) const { return 0.; }
  Real GetCvMassRef(int n) const { return 0.; }
  Real GetInvMu(int n) const { return 0.; }

 protected:
  void _rk4_integrate_lnp(Real dlnp, std::string method, Real adlnTdlnP) const;

  void _rk4_integrate_z(Real dlnp, std::string method, Real grav,
                        Real adlnTdlnP) const;

  template <typename T>
  Real _molar_entropy(T w, int i) const;

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
