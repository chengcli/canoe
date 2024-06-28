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
class Kinetics;
}  // namespace Cantera

using IndexPair = std::pair<int, int>;
using IndexSet = std::vector<int>;

using RealArray3 = std::array<Real, 3>;
using RealArrayX = std::vector<Real>;

enum { MAX_REACTANT = 3 };

using ReactionIndx = std::array<int, MAX_REACTANT>;
using ReactionStoi = std::array<int, MAX_REACTANT>;
using ReactionInfo = std::pair<ReactionIndx, ReactionStoi>;

void read_thermo_property(Real var[], char const name[], int len, Real v0,
                          ParameterInput *pin);
Real saha_ionization_electron_density(Real T, Real num, Real ion_ev);

//! Ideal saturation vapor pressure
//! $p^* = p^r\exp[\beta(1-1/t)-\delta\lnt]$
//! $p^r$ is the reference pressure, usually choosen to be the triple point
//! pressure \n $t=T/T^r$ is the dimensionless temperature. $T^r$ is the
//! reference temperature. \n Similar to $p^r$, $T^r$ is usually choosen to be
//! the triple point temperature
//! \return $p^*$ [pa]
inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real delta) {
  return p3 * exp((1. - 1. / t) * beta - delta * log(t));
}

// Thermodynamic variables are ordered in the array as the following
// 0: dry air (non-condensible gas)
// 1..NVAPOR: moist air (condensible gas)
// NVAPOR+1..NVAPOR+NCLOUD: clouds

// 0 is a special buffer place for cloud in equilibrium with vapor at the same
// temperature

class Thermodynamics {
 protected:
  //! Constructor for class sets up the initial conditions
  //! Protected ctor access thru static member function Instance
  Thermodynamics() {}
  static Thermodynamics *fromLegacyInput(ParameterInput *pin);
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

  //! Ideal gas constant of dry air in [J/(kg K)]
  //! \return $R_d=\hat{R}/\mu_d$
  Real GetRd() const { return Rd_; }

  //! reference adiabatic index of dry air [1]
  Real GetGammadRef() const { return gammad_ref_; }

  //! Molecular weight [kg/mol]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\mu$
  Real GetMu(int n) const { return mu_[n]; }

  //! Inverse molecular weight [mol/kg]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\mu$
  Real GetInvMu(int n) const { return inv_mu_[n]; }

  //! Ratio of molecular weights [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\epsilon_i=\mu_i/\mu_d$
  Real GetMuRatio(int n) const { return mu_ratio_[n]; }

  //! const pointer to mu_ratio_
  Real const *GetMuRatio() const { return &mu_ratio_[0]; }

  //! Ratio of molecular weights [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\epsilon_i^{-1} =\mu_d/\mu_i$
  Real GetInvMuRatio(int n) const { return inv_mu_ratio_[n]; }

  //! Ratio of specific heat capacity [J/(kg K)] at constant volume [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $c_{v,i}/c_{v,d}$
  Real GetCvRatioMass(int n) const { return cv_ratio_mass_[n]; }

  //! const pointer to cv_ratio_mass_
  Real const *GetCvRatioMass() const { return cv_ratio_mass_.data(); }

  //! Ratio of specific heat capacity [J/(mol K)] at constant volume [1]
  //! \param[in] n the index of the thermodynamic species
  //! \return $\hat{c}_{v,i}/\hat{c}_{v,d}$
  Real GetCvRatioMole(int n) const { return cv_ratio_mole_[n]; }

  //! const pointer to cv_ratio_mole_
  Real const *GetCvRatioMole() const { return cv_ratio_mole_.data(); }

  std::shared_ptr<Cantera::Kinetics> Kinetics() const { return kinetics_; }

  size_t GetSpeciesId(std::string const &name) const;

  //! Reference specific heat capacity [J/(kg K)] at constant volume
  Real GetCvMassRef(int n) const {
    Real cvd = Rd_ / (gammad_ref_ - 1.);
    return cv_ratio_mass_[n] * cvd;
  }

  //! Ratio of specific heat capacity [J/(kg K)] at constant pressure
  //! \return $c_{p,i}/c_{p,d}$
  Real GetCpRatioMass(int n) const { return cp_ratio_mass_[n]; }

  //! const pointer to cp_ratio_mass_
  Real const *GetCpRatioMass() const { return cp_ratio_mass_.data(); }

  //! Ratio of specific heat capacity [J/(mol K)] at constant pressure
  //! \return $\hat{c}_{p,i}/\hat{c}_{p,d}$
  Real GetCpRatioMole(int n) const { return cp_ratio_mole_[n]; }

  //! const pointer to cp_ratio_mole_
  Real const *GetCpRatioMole() const { return cp_ratio_mole_.data(); }

  Real GetCpMassRef(int n) const {
    Real cpd = Rd_ * gammad_ref_ / (gammad_ref_ - 1.);
    return cp_ratio_mass_[n] * cpd;
  }

  //! Thermodnamic equilibrium at current TP
  //! \param[in,out] qfrac mole fraction representation of air parcel
  void EquilibrateTP() const;

  //! Thermodnamic equilibrium at current UV
  void EquilibrateUV() const;

  void EquilibrateSP(double P) const;

  void FromPrimitive(MeshBlock const *pmb, int k, int j, int i) const;
  void ToPrimitive(MeshBlock *pmb, int k, int j, int i) const;

  void FromConserved(MeshBlock const *pmb, int k, int j, int i) const;
  void ToConserved(MeshBlock *pmb, int k, int j, int i) const;

  //! \brief Calculate potential temperature from primitive variable
  //!
  //! $\theta = T(\frac{p_0}{p})^{\chi}$
  //! \return $\theta$
  Real PotentialTemp(MeshBlock *pmb, Real p0, int k, int j, int i) const {
    auto &w = pmb->phydro->w;
    return GetTemp(pmb, k, j, i) *
           pow(p0 / w(IPR, k, j, i), GetChi(pmb, k, j, i));
  }

  //! \brief Calculate equivalent potential temperature from primitive variable
  //!
  //! $\theta_e = T(\frac{p}{p_d})^{Rd/(cpd + cl r_t} \exp(\frac{L_v q_v}{c_p
  //! T})$
  Real EquivalentPotentialTemp(MeshBlock *pmb, Real p0, int v, int k, int j,
                               int i) const;

  //! \brief Calculate temperature from primitive variable
  //!
  //! $T = p/(\rho R) = p/(\rho \frac{R}{R_d} Rd)$
  //! \return $T$
  Real GetTemp(MeshBlock const *pmb, int k, int j, int i) const {
    auto &w = pmb->phydro->w;
    return w(IPR, k, j, i) / (w(IDN, k, j, i) * Rd_ * RovRd(pmb, k, j, i));
  }

  //! \brief Mean molecular weight
  //!
  //! $mu = \mu_d (1 + \sum_i q_i (\epsilon_i - 1))$
  //! \return $mu$
  Real GetMu(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief Inverse of the mean molecular weight (no cloud)
  //!
  //! $ \frac{R}{R_d} = \frac{\mu_d}{\mu}$
  //! \return $1/\mu$
  //! Eq.16 in Li2019
  Real RovRd(MeshBlock const *pmb, int k, int j, int i) const {
    Real feps = 1.;
    auto &w = pmb->phydro->w;
#pragma omp simd reduction(+ : feps)
    for (int n = 1; n <= NVAPOR; ++n)
      feps += w(n, k, j, i) * (inv_mu_ratio_[n] - 1.);
    return feps;
  }

  //! \brief Specific heat capacity [J/(kg K)] of the air parcel at constant
  //! volume
  //!
  //! $c_v = c_{v,d}*(1 + \sum_i (q_i*(\hat{c}_{v,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{v,i}
  //!   - 1.)))$
  //! \return $c_v$
  Real GetCvMass(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief Specific heat capacity [J/(kg K)] of the air parcel at constant
  //! pressure
  //!
  //! $c_p = c_{p,d}*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{p,i}
  //!   - 1.)))$
  //! \return $c_p$
  Real GetCpMass(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief Adiabatic index
  //!
  //! $\chi = \frac{R}{c_p}$
  //! \return $\chi$
  Real GetChi(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief Polytropic index
  //!
  //! $\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  Real GetGamma(MeshBlock const *pmb, int k, int j, int i) const;

  //! Pressure from conserved variables
  //! \return $p$
  Real GetPres(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief specific enthalpy [J/kg] of the air parcel
  //!
  //! $h = c_{p,d}*T*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))$
  //! $  = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi}
  //! - 1.)))$
  Real GetEnthalpyMass(MeshBlock const *pmb, int k, int j, int i) const;

  //! \brief Moist static energy
  //!
  //! $h_s = c_{pd}T + gz + L_vq_v + L_s\sum_i q_i$
  //! \return $h_s$
  Real MoistStaticEnergy(MeshBlock const *pmb, Real gz, int k, int j,
                         int i) const;

  //! \brief Relative humidity
  Real RelativeHumidity(MeshBlock const *pmb, int n, int k, int j, int i) const;

 protected:
  //! custom functions for vapor (to be overridden by user mods)
  void enrollVaporFunctions();

 protected:
  std::shared_ptr<Cantera::Kinetics> kinetics_;

  //! ideal gas constant of dry air in J/kg
  Real Rd_;

  //! reference polytropic index of dry air
  Real gammad_ref_;

  //! ratio of mean molecular weights
  std::array<Real, Size> mu_ratio_;

  //! inverse ratio of mean molecular weights
  std::array<Real, Size> inv_mu_ratio_;

  //! molecular weight [kg/mol]
  std::array<Real, Size> mu_;

  //! inverse molecular weight [mol/kg]
  std::array<Real, Size> inv_mu_;

  //! ratio of specific heat capacities [J/mol] at constant pressure
  std::array<Real, Size> cp_ratio_mole_;

  //! ratio of specific heat capacities [J/mol] at constant pressure
  std::array<Real, Size> cp_ratio_mass_;

  //! ratio of specific heat capacities [J/mol] at constant volume
  std::array<Real, Size> cv_ratio_mole_;

  //! ratio of specific heat capacities [J/mol] at constant volume
  std::array<Real, Size> cv_ratio_mass_;

  //! pointer to the single Thermodynamics instance
  static Thermodynamics *mythermo_;

 private:
  // scratch space
  mutable std::array<Real, Size> mixr_;
};

#endif  // SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
