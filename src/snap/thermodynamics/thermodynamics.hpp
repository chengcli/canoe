#ifndef SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
#define SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_

// C/C++
#include <array>
#include <cfloat>
#include <iosfwd>
#include <memory>
#include <set>

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <constants.hpp>
#include <variable.hpp>

class MeshBlock;
class ParameterInput;

using CloudIndexSet = std::vector<int>;
using SatVaporPresFunc = Real (*)(Variable const &, int i, int j);

Real NullSatVaporPres(Variable const &, int, int);

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

class Thermodynamics {
 protected:
  enum { NPHASE = 3 };

  //! Constructor for class sets up the initial conditions
  //! Protected ctor access thru static member function Instance
  Thermodynamics() {}
  explicit Thermodynamics(YAML::Node &node);

 public:
  enum { Size = 1 + NVAPOR + NCLOUD };

  enum class Method {
    ReversibleAdiabat = 0,
    PseudoAdiabat = 1,
    DryAdiabat = 2,
    Isothermal = 3,
    NeutralStability = 4
  };

  // member functions
  ~Thermodynamics();

  //! Return a pointer to the one and only instance of Thermodynamics
  static Thermodynamics const *GetInstance();

  static Thermodynamics const *InitFromAthenaInput(ParameterInput *pin);

  //! Destroy the one and only instance of Thermodynamics
  static void Destroy();

  //! Ideal gas constant of dry air in [J/(kg K)]
  //! \return $R_d=\hat{R}/\mu_d$
  Real GetRd() const { return Rd_; }

  Real GetGammad(Variable const &var) const;

  Real GetGammadRef() const { return gammad_ref_; }

  //! Molecular weight
  //! \return $\mu$
  Real GetMu(int n) const { return mu_[n]; }

  //! Inverse molecular weight
  //! \return $\mu$
  Real GetInvMu(int n) const { return inv_mu_[n]; }

  //! Ratio of molecular weights
  //! \return $\epsilon_i=\mu_i/\mu_d$
  Real GetMuRatio(int n) const { return mu_ratio_[n]; }

  //! Ratio of molecular weights
  //! \return $\epsilon_i^{-1} =\mu_d/\mu_i$
  Real GetInvMuRatio(int n) const { return inv_mu_ratio_[n]; }

  //! Ratio of specific heat capacity [J/(kg K)] at constant volume
  //! \return $c_{v,i}/c_{v,d}$
  Real GetCvRatioMass(int n) const { return cv_ratio_mass_[n]; }

  //! Ratio of specific heat capacity [J/(mol K)] at constant volume
  //! \return $\hat{c}_{v,i}/\hat{c}_{v,d}$
  Real GetCvRatioMole(int n) const { return cv_ratio_mole_[n]; }

  int GetCloudIndex(int i, int j) const { return cloud_index_set_[i][j]; }

  //! Specific heat capacity [J/(kg K)] at constant volume
  //! $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
  //! $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
  //! \return $c_{v,i}$
  Real GetCvMass(Variable const &qfrac, int n) const {
    Real cvd = Rd_ / (GetGammad(qfrac) - 1.);
    return cv_ratio_mass_[n] * cvd;
  }

  Real GetCvMassRef(int n) const {
    Real cvd = Rd_ / (gammad_ref_ - 1.);
    return cv_ratio_mass_[n] * cvd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant volume
  //! \return $\hat{c}_v$
  Real GetCvMole(Variable const &qfrac, int n) const {
    return GetCvMass(qfrac, n) * mu_[n];
  }

  //! Ratio of specific heat capacity [J/(kg K)] at constant pressure
  //! \return $c_{p,i}/c_{p,d}$
  Real GetCpRatioMass(int n) const { return cp_ratio_mass_[n]; }

  //! Ratio of specific heat capacity [J/(mol K)] at constant pressure
  //! \return $\hat{c}_{p,i}/\hat{c}_{p,d}$
  Real GetCpRatioMole(int n) const { return cp_ratio_mole_[n]; }

  //! Specific heat capacity [J/(kg K)] at constant pressure
  //! $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
  //! $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
  //! \return $c_p$
  Real GetCpMass(Variable const &qfrac, int n) const {
    Real gammad = GetGammad(qfrac);
    Real cpd = Rd_ * gammad / (gammad - 1.);
    return cp_ratio_mass_[n] * cpd;
  }

  Real GetCpMassRef(int n) const {
    Real cpd = Rd_ * gammad_ref_ / (gammad_ref_ - 1.);
    return cp_ratio_mass_[n] * cpd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant pressure
  //! \return $\hat{c}_v$
  Real GetCpMole(Variable const &qfrac, int n) const {
    return GetCpMass(qfrac, n) * mu_[n];
  }

  //! Adiabatic index
  //! \return $\chi = \frac{R}{c_p}$
  Real GetChi(Variable const &qfrac) const;

  //! Temperature dependent specific latent energy [J/kg] of condensates at
  //! constant volume $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T - T^r)$
  //! $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMass(int n, Real temp) const {
    return latent_energy_mass_[n] - delta_[n] * Rd_ * inv_mu_ratio_[n] * temp;
  }

  Real GetLatentEnergyMass(int n) const { return latent_energy_mass_[n]; }

  //! Temperature dependent specific latent energy [J/mol] of condensates at
  //! constant volume
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMole(int n, Real temp) const {
    return GetLatentEnergyMass(n, temp) * mu_[n];
  }

  Real GetLatentEnergyMole(int n) const { return latent_energy_mole_[n]; }

  Real GetLatentHeatMole(int i, std::vector<Real> const &rates,
                         Real temp) const;

  Real GetLatentHeatMass(int i, std::vector<Real> const &rates,
                         Real temp) const {
    return GetLatentHeatMole(i, rates, temp) * inv_mu_[i];
  }

  //! Calculate the equilibrium mole transfer between vapor and cloud
  //! \param L_ov_cv L/cv evaluated at current temperature
  //! \return molar fraction change of vapor to cloud
  std::vector<Real> TryEquilibriumTP(Variable const &qfrac, int ivapor,
                                     Real cv_hat = 0.,
                                     bool misty = false) const;

  //! Construct an 1d atmosphere
  //! @param method choose from [reversible, pseudo, dry, isothermal]
  void Extrapolate(Variable *qfrac, Real dzORdlnp, Method method,
                   Real grav = 0., Real userp = 0.) const;

  //! Adjust to the maximum saturation state conserving internal energy
  void SaturationAdjustment(Variable *qfrac) const;

  //! Potential temperature
  //!$\theta = T(\frac{p_0}{p})^{\chi}$
  //! \return $\theta$
  Real PotentialTemp(MeshBlock *pmb, Real p0, int k, int j, int i) const {
    auto &w = pmb->phydro->w;
    return GetTemp(pmb, k, j, i) *
           pow(p0 / w(IPR, k, j, i), GetChi(pmb, k, j, i));
  }

  //! Equivalent potential temperature
  //! $\theta_e = T(\frac{p}{p_d})^{Rd/(cpd + cl r_t} \exp(\frac{L_v q_v}{c_p
  //! T})$
  Real EquivalentPotentialTemp(MeshBlock *pmb, Real p0, int v, int k, int j,
                               int i) const;

  //! Temperature
  //!$T = p/(\rho R) = p/(\rho \frac{R}{R_d} Rd)
  //! \return $T$
  Real GetTemp(MeshBlock *pmb, int k, int j, int i) const {
    auto &w = pmb->phydro->w;
    return w(IPR, k, j, i) / (w(IDN, k, j, i) * Rd_ * RovRd(pmb, k, j, i));
  }

  //! Mean molecular weight
  //! $mu = \mu_d (1 + \sum_i q_i (\epsilon_i - 1))$
  //! \return $mu$
  Real GetMu(MeshBlock *pmb, int k, int j, int i) const;

  //! Inverse of the mean molecular weight
  //! $ \frac{R}{R_d} = \frac{\mu_d}{\mu}$
  //! \return $1/\mu$
  Real RovRd(MeshBlock *pmb, int k, int j, int i) const;

  Real RovRd(Variable const &qfrac) const;

  //! Specific heat capacity [J/(kg K)] of the air parcel at constant volume
  //! c_v = c_{v,d}*(1 + \sum_i (q_i*(\hat{c}_{v,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{v,i} - 1.)))
  //! \return $c_v$
  Real GetCvMass(MeshBlock *pmb, int k, int j, int i) const;

  //! Specific heat capacity [J/(kg K)] of the air parcel at constant pressure
  //! c_p = c_{p,d}*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))
  //! \return $c_p$
  Real GetCpMass(MeshBlock *pmb, int k, int j, int i) const;

  //! Adiabatic index
  //!$\chi = \frac{R}{c_p}$
  //! \return $\chi$
  Real GetChi(MeshBlock *pmb, int k, int j, int i) const;

  //! Polytropic index
  //!$\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  Real GetGamma(MeshBlock *pmb, int k, int j, int i) const;

  //! Pressure from conserved variables
  //! \return $p$
  Real GetPres(MeshBlock *pmb, int k, int j, int i) const;

  //! Enthalpy
  //!$h = c_{p,d}*T*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))$
  //!$  = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))$
  Real GetEnthalpyMass(MeshBlock *pmb, int k, int j, int i) const;

  //! Moist static energy
  //!$h_s = c_{pd}T + gz + L_vq_v + L_s\sum_i q_i$
  //! \return $h_s$
  Real MoistStaticEnergy(MeshBlock *pmb, Real gz, int k, int j, int i) const;

  //! Relative humidity
  //! $H_i = \frac{e_i}{e_i^s}$
  //! \return $H$
  Real RelativeHumidity(MeshBlock *pmb, int n, int k, int j, int i) const;

  Real RelativeHumidity(Variable const &qfrac, int n) const {
    auto rates = TryEquilibriumTP(qfrac, n, 0., true);
    return qfrac.w[n] / (qfrac.w[n] + rates[0]);
  }

 protected:
  //! update T/P
  void updateTPConservingU(Variable *qfrac, Real rmole, Real umole) const;

  Real getInternalEnergyMole(Variable const &qfrac) const;

  Real getDensityMole(Variable const &qfrac) const {
    Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
    for (int n = 0; n < NCLOUD; ++n) qgas += -qfrac.c[n];
    return qfrac.w[IPR] / (Constants::Rgas * qfrac.w[IDN] * qgas);
  }

  void setTotalEquivalentVapor(Variable *qfrac, int i) const {
    for (auto &j : cloud_index_set_[i]) {
      qfrac->w[i] += qfrac->c[j];
      qfrac->c[j] = 0.;
    }
  }

  //! Calculate moist adiabatic temperature gradient
  //! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
  //! \return $\Gamma_m$
  Real calDlnTDlnP(Variable const &qfrac, Real latent[]) const;

  void rk4IntegrateLnp(Variable *qfrac, Real dlnp, Method method,
                       Real adlnTdlnP) const;
  void rk4IntegrateZ(Variable *qfrac, Real dlnp, Method method, Real grav,
                     Real adlnTdlnP) const;

  void enrollVaporFunctions(ParameterInput *pin);
  void enrollVaporFunctionsEarth();
  void enrollVaporFunctionsGiants();
  void enrollVaporFunctionsJupiterJuno();

 private:
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

  //! latent energy at constant volume [J/kg]
  //! $L_i^r+\Delta c_{ij}T^r == \beta_{ij}\frac{R_d}{\epsilon_i}T^r$
  std::array<Real, Size> latent_energy_mass_;

  //! dimensionless latent heat
  //! $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$
  //! $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$
  //! and the condensate $j$ $R_i=\hat{R}/m_i=R_d/\epsilon_i$ $T^r$ is the
  //! triple point temperature
  //! $\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r} $
  std::array<Real, Size> beta_;

  //! dimensionless differences in specific heat capacity
  //! $\delta_{ij} = \frac{c_{ij} - c_{p,i}}{R_i}$
  //! $c_{ij} - c_{p,j}$ is the difference in specific heat
  //! capacity at constant pressure. $c_{ij}$ is the heat capacity of
  //! condensate $j$ of vapor $i$
  std::array<Real, Size> delta_;

  //! triple point temperature [K]
  std::array<Real, 1 + NVAPOR> t3_;

  //! triple point pressure [pa]
  std::array<Real, 1 + NVAPOR> p3_;

  //! saturation vapor pressure function
  std::vector<std::vector<SatVaporPresFunc>> svp_func_;

  //! cloud index set
  std::vector<CloudIndexSet> cloud_index_set_;

  //! Pointer to the single Application instance
  static Thermodynamics *mythermo_;

  //! other parameters
  int sa_max_iter_ = 5;
  Real sa_relax_ = 0.8;
  Real sa_ftol_ = 0.01;
};

#endif  // SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
