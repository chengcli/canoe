#ifndef SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
#define SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_

// C/C++
#include <array>
#include <cfloat>
#include <iosfwd>
#include <memory>

// canoe
#include <configure.hpp>
#include <constants.hpp>
#include <variable.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/mesh/mesh.hpp>

class MeshBlock;
class ParameterInput;

// dynamic variables are ordered in the array as the following
// 0: dry air (non-condensible)
// 1+iphase*NVAPOR:1+(iphase+1)*NVAPOR
// iphase = 0 - moist air (condensible)
// iphase = 1 - primary condensible species
// iphase = 2..(N-1) - all other condensible species

class Thermodynamics {
  friend std::ostream &operator<<(std::ostream &os, Thermodynamics const &my);

 public:
  enum { Size = 1 + NVAPOR + NCLOUD };

  enum class Method {
    ReversibleAdiabat,
    PseudoAdiabat,
    DryAdiabat,
    Isothermal
  };

  // member functions
  Thermodynamics(MeshBlock *pmb, ParameterInput *pin);
  ~Thermodynamics();

  void EnrollSatVaporPresFuncs() {
    if (iv == AMMONIA_VAPOR_ID)
      s = sat_vapor_p_NH3_BriggsS(q[IDN]);
    else if (iv == WATER_VAPOR_ID)
      s = sat_vapor_p_H2O_BriggsS(q[IDN]);
    else
      s = SatVaporPresIdeal(t, p3, beta, delta);
    s /= q[IPR];
  }

  //! Calculate the equilibrium mole transfer between vapor and cloud
  //! \param L_ov_cv L/cv evaluated at current temperature
  //! \return molar fraction change of vapor to cloud
  Real TryEquilibriumTP(Variable const &qfrac, int ivapor, Real L_ov_cv = 0.,
                        bool no_cloud = false) const;

  Real ExecuteEquilibriumTP(Variable const &qfrac, int ivapor);

  //! Calculate moist adiabatic temperature gradient
  //! $\Gamma_m = (\frac{d\ln T}{d\ln P})_m$
  //! \return $\Gamma_m$
  Real CalDlnTDlnP(Variable const &qfrac, Real gammad, int isat[]) const;

  //! Ideal gas constant of dry air in [J/(kg K)]
  //! \return $R_d=\hat{R}/\mu_d$
  Real GetRd() const { return Rd_; }

  //! Molecular weight
  //! \return $\mu$
  Real GetMu(int n) const { return mu_[n]; }

  //! Mean molecular weight
  //! $mu = \mu_d (1 + \sum_i q_i (\epsilon_i - 1))$
  //! \return $mu$
  Real GetMu(int k, int j, int i) const;

  //! Ratio of molecular weights
  //! \return $\epsilon_i=\mu_i/\mu_d$
  Real GetMassRatio(int n) const { return mu_ratio_[n]; }

  //! Inverse of the mean molecular weight
  //! $ \frac{R}{R_d} = \frac{\mu_d}{\mu}$
  //! \return $1/\mu$
  Real RovRd(int k, int j, int i) const;

  //! Ratio of specific heat capacity [J/(kg K)] at constant volume
  //! \return $c_{v,i}/c_{v,d}$
  Real GetCvRatioMass(int n) const { return cv_ratio_mass_[n]; }

  //! Ratio of specific heat capacity [J/(mol K)] at constant volume
  //! \return $\hat{c}_{v,i}/\hat{c}_{v,d}$
  Real GetCvRatioMole(int n) const { return cv_ratio_mole_[n]; }

  //! Specific heat capacity [J/(kg K)] at constant volume
  //! $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
  //! $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
  //! \return $c_{v,i}$
  Real GetCvMass(int n) const {
    Real cvd = Rd_ / (pmy_block_->peos->GetGamma() - 1.);
    return cv_ratio_mass_[n] * cvd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant volume
  //! \return $\hat{c}_v$
  Real GetCvMole(int n) const { return GetCvMass(n) * mu_[n]; }

  //! Specific heat capacity [J/(kg K)] of the air parcel at constant volume
  //! c_v = c_{v,d}*(1 + \sum_i (q_i*(\hat{c}_{v,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{v,i} - 1.)))
  //! \return $c_v$
  Real GetCvMass(int k, int j, int i) const;

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
  Real GetCpMass(int n) const {
    Real gamma = pmy_block_->peos->GetGamma();
    Real cpd = Rd_ * gamma / (gamma - 1.);
    return cp_ratio_mass_[n] * cpd;
  }

  //! Specific heat capacity [J/(mol K)] of the air parcel at constant pressure
  //! \return $\hat{c}_v$
  Real GetCpMole(int n) const { return GetCpMass(n) * mu_[n]; }

  //! Specific heat capacity [J/(kg K)] of the air parcel at constant pressure
  //! c_p = c_{p,d}*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{p,i} - 1.)))
  //! \return $c_p$
  Real GetCpMass(int k, int j, int i) const;

  //! Temperature dependent specific latent energy [J/kg] of condensates at
  //! constant volume $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T - T^r)$
  //! $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMass(int n, Real temp) const {
    return latent_[n] - delta_[n] * Rd_ / mu_ratio_[j] * temp;
  }

  //! Temperature dependent specific latent energy [J/mol] of condensates at
  //! constant volume
  //! \return $L_{ij}(T)$
  Real GetLatentEnergyMole(int n, Real temp) const {
    return GetLatentMass(n, temp) * mu_[n];
  }

  Real GetLatentHeat(int n, Real temp) const {
    // int nc = iv + NVAPOR;
    int nc = q[IDN] > t3[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
    Real latent = beta[nc] * t3[iv] / q[IDN] - delta[nc];
  }

  //! Dimensionless latent heat
  //! $\beta$ parameter of a condensate
  //! $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$
  //! $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$
  //! and the condensate $j$ $R_i=\hat{R}/m_i=R_d/\epsilon_i$ $T^r$ is the
  //! triple point temperature
  //! \return $\beta_{ij}$
  Real GetBeta(int n) const { return beta_[n]; }

  //! Dimensionless heat capacity difference
  //! $\delta$ parameter of a condensate
  //! $\delta_{ij} = \frac{\Delta c_{ij}}{R_i}$
  //! $\Delta c_{ij} = c_{ij} - c_{p,j}$ is the difference in specific heat
  //! capacity at constant temperature. $c_{ij}$ is the heat capacity of
  //! condensate $j$ of vapor $i$
  //! \return $\delta_{ij}$
  Real GetDelta(int n) const { return delta_[n]; }

  //! Adiabatic index
  //!$\chi = \frac{R}{c_p}$
  //! \return $\chi$
  Real GetChi(int k, int j, int i) const;

  //! Polytropic index
  //!$\gamma = \frac{c_p}{c_v}$
  //! \return $\gamma$
  Real GetGamma(int k, int j, in i) const;

  //! Potential temperature
  //!$\theta = T(\frac{p_0}{p})^{\chi}$
  //! \return $\theta$
  Real PotentialTemp(Real p0, int k, int j, int i) const {
    auto &w = pmy_block_->phydro->w;
    return GetTemp(k, j, i) * pow(p0 / w(IPR, k, j, i), GetChi(k, j, i));
  }

  //! Temperature
  //!$T = p/(\rho R) = p/(\rho \frac{R}{R_d} Rd)
  //! \return $T$
  Real GetTemp(int k, int j, int i) const {
    auto &w = pmy_block_->phydro->w;
    return w(IPR, k, j, i) / (w(IDN, k, j, i) * Rd_ * RovRd(k, j, i));
  }

  //! Pressure from conserved variables
  //! \return $p$
  Real GetPres(int k, int j, int i) const;

  //! Construct an 1d atmosphere
  //! @param method choose from [reversible, pseudo, dry, isothermal]
  void ConstructAtmosphere(Variable *prim, Real dzORdlnp, Adiabat method,
                           Real userp) const;

  //! Enthalpy
  //!$h = c_{pd}*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))$
  //!$  = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))$
  Real GetEnthalpyMass(int k, int j, int i) const;

  //! Moist static energy
  //!$h_s = c_{pd}T + gz + L_vq_v + L_s\sum_i q_i$
  //! \return $h_s$
  Real MoistStaticEnergy(Real gz, int k, int j, int i);

  //! Saturation surplus for vapors can be both positive and negative
  //! positive value represents supersaturation
  //! negative value represents saturation deficit
  void SaturationSurplus(Real dv[], Variable const &v,
                         Variable::Type vtype) const;

  //! Relative humidity
  //! $H = \frac{e_i}{e_i^s}$
  //! \return $H$
  Real RelativeHumidity(int n, int k, int j, int i) const;

  //! Retrieve saturation vapor pressure function
  SatVaporPresFunc GetSatVaporPresFunc(int ivapor, int cloud) const;

 private:
  MeshBlock *pmy_block_;

  //! ideal gas constant of dry air in J/kg
  Real Rd_;

  //! ratio of mean molecular weights
  std::array<Real, Size> mu_ratio_;

  //! molecular weight [kg/mol]
  std::array<Real, Size> mu_;

  //! ratio of specific heat capacities [J/mol] at constant pressure
  std::array<Real, Size> cp_ratio_mole_;

  //! ratio of specific heat capacities [J/mol] at constant pressure
  std::array<Real, Size> cp_ratio_mass_;

  //! ratio of specific heat capacities [J/mol] at constant volume
  std::array<Real, Size> cv_ratio_mole_;

  //! ratio of specific heat capacities [J/mol] at constant volume
  std::array<Real, Size> cv_ratio_mass_;

  //! latent heat at constant pressure [J/mol]
  std::array<Real, Size> latent_heat_mole_;

  //! latent heat at constant pressure [J/kg]
  //! $L_i^r+\Delta c_{ij}T^r == \beta_{ij}\frac{R_d}{\epsilon_i}T^r$
  std::array<Real, Size> latent_heat_mass_;

  //! dimensionless latent heat
  //! $\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r} $
  std::array<Real, Size> beta_;

  //! dimensionless differences in specific heat capacity
  //! $(c_{ij} - c_{p,i})/R_i$
  std::array<Real, Size> delta_;

  //! triple point temperature [K]
  std::array<Real, 1 + NVAPOR> t3_;

  //! triple point pressure [pa]
  std::array<Real, 1 + NVAPOR> p3_;

  //! saturation vapor pressure function
  std::vector<SatVaporPresFunc> svp_func_;

 protected:
  void updateGammad(Variable const &var);

  void rk4IntegrateLnp(Real q[], int isat[], Real const rcp[],

  Real setTotalEquivalentVapor(qfrac) {
    q[iv] += q[NHYDRO + iv - 1] + q[NHYDRO + NVAPOR + iv - 1];
    q[NHYDRO + iv - 1] = 0.;
    q[NHYDRO + NVAPOR + iv - 1] = 0.;
  }
};

using ThermodynamicsPtr = std::shared_ptr<Thermodynamics>;

#endif  // SRC_SNAP_THERMODYNAMICS_THERMODYNAMICS_HPP_
