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

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>

class MeshBlock;
class ParameterInput;

using IndexPair = std::pair<int, int>;
using IndexSet = std::vector<int>;

using RealArray3 = std::array<Real, 3>;
using RealArrayX = std::vector<Real>;

using SatVaporPresFunc1 = Real (*)(AirParcel const &, int i, int j);
using SatVaporPresFunc2 = Real (*)(AirParcel const &, int i, int j, int k);

//! \todo(CLI): move to configure.hpp
enum { MAX_REACTANT = 3 };
using ReactionIndx = std::array<int, MAX_REACTANT>;
using ReactionStoi = std::array<int, MAX_REACTANT>;
using ReactionInfo = std::pair<ReactionIndx, ReactionStoi>;

Real NullSatVaporPres1(AirParcel const &, int, int);
Real NullSatVaporPres2(AirParcel const &, int, int, int);

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
  enum { NPHASE = NPHASE_LEGACY };

  //! Constructor for class sets up the initial conditions
  //! Protected ctor access thru static member function Instance
  Thermodynamics() {}
  static Thermodynamics *fromLegacyInput(ParameterInput *pin);
  static Thermodynamics *fromYAMLInput(YAML::Node const &node);

 public:
  using SVPFunc1Container = std::vector<std::vector<SatVaporPresFunc1>>;
  using SVPFunc2Container = std::map<IndexPair, SatVaporPresFunc2>;

  enum { Size = 1 + NVAPOR + NCLOUD };

  static constexpr Real RefTemp = 300.;
  static constexpr Real RefPres = 1.e5;

  //! thermodynamics input key in the input file [thermodynamics_config]
  static const std::string input_key;

  // member functions
  ~Thermodynamics();

  //! Return a pointer to the one and only instance of Thermodynamics
  static Thermodynamics const *GetInstance();

  static Thermodynamics const *InitFromYAMLInput(YAML::Node const &node);

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

  //! \return the index of the cloud
  //! \param[in] i the index of the vapor
  //! \param[in] j the sequential index of the cloud
  int GetCloudIndex(int i, int j) const { return cloud_index_set_[i][j]; }

  //! const pointer to cloud_index_set_
  IndexSet const *GetCloudIndexSet() const { return cloud_index_set_.data(); }

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

  //! \brief Temperature dependent specific latent energy [J/kg] of condensates
  //! at constant volume $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T -
  //! T^r)$
  //!
  //! $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMass(int n, Real temp = 0.) const {
    return latent_energy_mass_[n] - delta_[n] * Rd_ * inv_mu_ratio_[n] * temp;
  }

  //! \brief Temperature dependent specific latent energy [J/mol] of condensates
  //! at constant volume
  //!
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMole(int n, Real temp = 0.) const {
    return GetLatentEnergyMass(n, temp) * mu_[n];
  }

  //! const pointer to latent_energy_mole_
  Real const *GetLatentEnergyMole() const { return &latent_energy_mole_[0]; }

  Real GetLatentHeatMole(int i, std::vector<Real> const &rates,
                         Real temp) const;

  Real GetLatentHeatMass(int i, std::vector<Real> const &rates,
                         Real temp) const {
    return GetLatentHeatMole(i, rates, temp) * inv_mu_[i];
  }

  //! \brief Calculate the equilibrium mole transfer by cloud reaction
  //! vapor -> cloud
  //!
  //! \param[in] qfrac mole fraction representation of air parcel
  //! \param[in] ivapor the index of the vapor
  //! \param[in] cv_hat $cv_hat$ molar heat capacity
  //! \param[in] misty if true, there is an infinite supple of cloud
  //! \return molar fraction change of vapor to cloud
  RealArrayX TryEquilibriumTP_VaporCloud(AirParcel const &qfrac, int ivapor,
                                         Real cv_hat = 0.,
                                         bool misty = false) const;

  //! \brief Calculate the equilibrium mole transfer by cloud reaction
  //! vapor1 + vapor2 -> cloud
  //!
  //! \param[in] air mole fraction representation of air parcel
  //! \param[in] ij the index pair of the vapor1 and vapor2
  //! \param[in] cv_hat $\har{cv}$ molar heat capacity
  //! \param[in] misty if true, there is an infinite supple of cloud
  //! \return molar fraction change of vapor1, vapor2 and cloud
  RealArray3 TryEquilibriumTP_VaporVaporCloud(AirParcel const &air,
                                              IndexPair ij, Real cv_hat = 0.,
                                              bool misty = false) const;

  //! Construct an 1d atmosphere
  //! \param[in,out] qfrac mole fraction representation of air parcel
  //! \param[in] dzORdlnp vertical grid spacing
  //! \param[in] method choose from [reversible, pseudo, dry, isothermal]
  //! \param[in] grav gravitational acceleration
  //! \param[in] userp user parameter to adjust the temperature gradient
  void Extrapolate(AirParcel *qfrac, Real dzORdlnp, std::string method,
                   Real grav = 0., Real userp = 0.) const;

  //! Thermodnamic equilibrium at current TP
  //! \param[in,out] qfrac mole fraction representation of air parcel
  void EquilibrateTP(AirParcel *qfrac) const;

  //! Adjust to the maximum saturation state conserving internal energy
  //! \param[in,out] ac mole fraction representation of a collection of air
  //! parcels
  void SaturationAdjustment(AirColumn &ac) const;

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

  //! latent energy at constant volume [J/mol]
  std::array<Real, Size> latent_energy_mole_;

  //! \brief latent energy at constant volume [J/kg]
  //!
  //! $L_i^r+\Delta c_{ij}T^r == \beta_{ij}\frac{R_d}{\epsilon_i}T^r$
  std::array<Real, Size> latent_energy_mass_;

  //! \brief Dimensionless latent heat
  //!
  //! $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$
  //! $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$
  //! and the condensate $j$ $R_i=\hat{R}/m_i=R_d/\epsilon_i$ $T^r$ is the
  //! triple point temperature
  //! $\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r}$
  std::array<Real, Size> beta_;

  //! \brief Dimensionless differences in specific heat capacity
  //!
  //! $\delta_{ij} = \frac{c_{ij} - c_{p,i}}{R_i}$
  //! $c_{ij} - c_{p,j}$ is the difference in specific heat
  //! capacity at constant pressure. $c_{ij}$ is the heat capacity of
  //! condensate $j$ of vapor $i$
  std::array<Real, Size> delta_;

  //! triple point temperature [K]
  std::array<Real, 1 + NVAPOR> t3_;

  //! triple point pressure [pa]
  std::array<Real, 1 + NVAPOR> p3_;

  //! saturation vapor pressure function: Vapor -> Cloud
  SVPFunc1Container svp_func1_;

  //! cloud index set
  std::vector<IndexSet> cloud_index_set_;

  //! saturation vapor pressure function: Vapor + Vapor -> Cloud
  SVPFunc2Container svp_func2_;

  //! reaction information map
  std::map<IndexPair, ReactionInfo> cloud_reaction_map_;

  //! pointer to the single Thermodynamics instance
  static Thermodynamics *mythermo_;

  //! maximum saturation adjustment iterations
  int sa_max_iter_ = 5;

  //! saturation adjustment relaxation parameter
  Real sa_relax_ = 0.8;

  //! saturation adjustment temperature tolerance
  Real sa_ftol_ = 0.01;
};

#endif  // SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
