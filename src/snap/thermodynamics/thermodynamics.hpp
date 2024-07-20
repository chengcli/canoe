#ifndef SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
#define SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_

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

class MeshBlock;
class ParameterInput;

namespace Cantera {
class ThermoPhase;
class Condensation;
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

  void UpdateThermoProperties(ThermoPhase &thermo);

  //! Ideal gas constant of dry air in [J/(kg K)]
  //! \return $R_d=\hat{R}/\mu_d$
  Real GetRd() const { return Rd_; }

  //! reference adiabatic index of dry air [1]
  Real GetGammad() const { return gammad_; }

  //! Construct an 1d atmosphere
  //! \param[in,out] qfrac mole fraction representation of air parcel
  //! \param[in] dzORdlnp vertical grid spacing
  //! \param[in] method choose from [reversible, pseudo, dry, isothermal]
  //! \param[in] grav gravitational acceleration
  //! \param[in] userp user parameter to adjust the temperature gradient
  template <typename T>
  void Extrapolate_inplace(T w, Real dzORdlnp, std::string method,
                           Real grav = 0., Real userp = 0.) const;

  //! Thermodnamic equilibrium at current TP
  //! \param[in,out] qfrac mole fraction representation of air parcel
  void EquilibrateTP() const;

  //! Thermodnamic equilibrium at current UV
  void EquilibrateUV() const;

  template <typename T>
  void SetStateFromPrimitive(T w) const;

  template <typename T>
  void SetStateFromConserved(T u) const;

 public:
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
    auto &w = pmb->phydro->w;
    return w[IPR] / (w[IDN] * Rd_ * RovRd(w));
  }

  //! Pressure from conserved variables
  //! \return $p$
  template <typename T>
  Real GetPres(T u, std::vector<Real> cos_theta) const;

  //! \brief Calculate potential temperature from primitive variable
  //!
  //! $\theta = T(\frac{p_0}{p})^{\chi}$
  //! \return $\theta$
  template <typename T>
  Real PotentialTemp(T w) const {
    return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
  }

  //! \brief Effective polytropic index
  //!
  //! Eq.63 in Li2019
  //! $\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  template <typename T>
  Real GetGamma(T w) const;

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

 protected:
  void _rk4_integrate_lnp(Real dlnp, std::string method, Real adlnTdlnP);

  void _rk4_integrate_z(Real dlnp, std::string method, Real grav,
                        Real adlnTdlnP);

 protected:
  std::shared_ptr<Cantera::Condensation> kinetics_;

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
//! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
//! \return $\Gamma_m$
Real cal_dlnT_dlnP(AirParcel const &qfrac, Real const *cp_ratio_mole,
                   Real const *latent);

#endif  // SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
