/** @file thermodynamics.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 25, 2021 18:04:45 UTC
 * @bug No known bugs.
 */

#ifndef THERMODYNAMICS_HPP
#define THERMODYNAMICS_HPP

// C/C++ headers
#include <cfloat>
#include <iosfwd>

// Athena headers
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/mesh/mesh.hpp>
#include <configure.hpp>

#include "moist_adiabat_funcs.hpp"

class MeshBlock;
class ParameterInput;

// dynamic variables are ordered in the array as the following
// 0: dry air (non-condensible)
// 1+iphase*NVAPOR:1+(iphase+1)*NVAPOR
// iphase = 0 - moist air (condensible)
// iphase = 1 - primary condensible species
// iphase = 2..(N-1) - all other condensible species

enum class Adiabat { reversible = 0, pseudo = 1, dry = 2, isothermal = 3 };

class Thermodynamics {
  friend std::ostream &operator<<(std::ostream &os, Thermodynamics const &my);

 public:
  // members
  static Real const Rgas;
  static Real const kBoltz;
  static Real const kBoltz_cgs;
  MeshBlock *pmy_block; /**< pointer to MeshBlock */

  // member functions
  Thermodynamics(MeshBlock *pmb, ParameterInput *pin);
  ~Thermodynamics() {}

  /*! Ideal gas constant of dry air in [J/(kg K)]
   * @return $R_d=\hat{R}/m_d$
   */
  Real GetRd() const { return Rd_; }

  /*! Ratio of specific heat capacity at constant volume
   * @param i index of the vapor
   * @return $c_{v,i}/c_{v,d}$
   */
  Real GetCvRatio(int i) const { return cv_ratios_[i]; }

  /*! Specific heat capacity at constant volume
   *
   * $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
   * $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
   * @param i index of the vapor
   * @return $c_v$ [J/(kg K)]
   */
  Real GetCv(int i) const {
    Real cvd = Rd_ / (pmy_block->peos->GetGamma() - 1.);
    return cv_ratios_[i] * cvd;
  }

  /*! Ratio of specific heat capacity at constant pressure
   * @param i index of the vapor
   * @return $c_{p,i}/c_{p,d}$
   */
  Real GetCpRatio(int i) const { return cp_ratios_[i]; }

  /*! Specific heat capacity at constant pressure
   *
   * $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
   * $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
   * @param i index of the vapor
   * @return $c_p$ [J/(kg K)]
   * */
  Real GetCp(int i) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real cpd = Rd_ * gamma / (gamma - 1.);
    return cp_ratios_[i] * cpd;
  }

  /*! Temperature dependent specific latent heat of condensates at constant
   * volume
   *
   * $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T - T^r)$
   * $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
   * @param j index of condensate
   * @param temp temperature (default to 0)
   * @return $L_{ij}(T)$ [J/kg]
   */
  Real GetLatent(int j, Real temp = 0.) const {
    return latent_[j] - delta_[j] * Rd_ / mu_ratios_[j] * temp;
  }

  /*! Ratio of molecular weights
   * @param i index of the vapor or the condensate
   * @return $\epsilon_i=m_i/m_d$
   */
  Real GetMassRatio(int i) const { return mu_ratios_[i]; }

  /*! $\beta$ parameter of a condensate
   *
   * $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$ \n
   * $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$
   * and the condensate $j$ \n
   * $R_i=\hat{R}/m_i=R_d/\epsilon_i$ \n
   * $T^r$ is the triple point temperature
   * @param j index of the condensate
   * @return $\beta_{ij}$
   */
  Real GetBeta(int j) const { return beta_[j]; }

  /*! $\delta$ parameter of a condensate
   *
   * $\delta_{ij} = \frac{\Delta c_{ij}}{R_i}$ \n
   * $\Delta c_{ij} = c_{ij} - c_{p,j}$ is the difference in specific heat
   * capacity at constant temperature.\n $c_{ij}$ is the heat capacity of
   * condensate $j$ of vapor $i$
   * @param j index of the condensate
   * @return $\delta_{ij}$
   */
  Real GetDelta(int j) const { return delta_[j]; }

  /*! Construct an 1d atmosphere
   * @param method choose from [reversible, pseudo, dry, isothermal]
   */
  void ConstructAtmosphere(Real **w, Real Ts, Real Ps, Real grav, Real dzORdlnp,
                           int len, Adiabat method, Real userp) const;

  // conversion functions
  //! Change mass mixing ratio to molar mixing ratio
  template <typename T1, typename T2>
  void PrimitiveToChemical(T1 c, T2 const w) const {
    // set molar mixing ratio
    Real sum = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      c[n] = w[n] / mu_ratios_[n];
      sum += w[n] * (1. / mu_ratios_[n] - 1.);
    }
    // set pressure, temperature, velocity
    c[IPR] = w[IPR];
    c[IDN] = w[IPR] / (w[IDN] * Rd_ * sum);
    c[IVX] = w[IVX];
    c[IVY] = w[IVY];
    c[IVZ] = w[IVZ];

    Real mols = c[IPR] / (c[IDN] * Rgas);
#pragma omp simd
    for (int n = 1; n <= NVAPOR; ++n) {
      c[n] *= mols / sum;
    }
  }

  //! Change molar mixing ratio to mass mixing ratio
  template <typename T1, typename T2>
  void ChemicalToPrimitive(T1 w, T2 const c) const {
    // set mass mixing ratio
    Real sum = 1., mols = c[IPR] / (c[IDN] * Rgas);
    for (int n = 1; n <= NVAPOR; ++n) {
      w[n] = c[n] / mols * mu_ratios_[n];
      sum += c[n] / mols * (mu_ratios_[n] - 1.);
    }
#pragma omp simd
    for (int n = 1; n <= NVAPOR; ++n) w[n] /= sum;

    // set pressure, density, velocity
    w[IPR] = c[IPR];
    w[IDN] = sum * c[IPR] / (c[IDN] * Rd_);
    w[IVX] = c[IVX];
    w[IVY] = c[IVY];
    w[IVZ] = c[IVZ];
  }

  //! Change density to molar mixing ratio
  template <typename T1, typename T2>
  void ConservedToChemical(T1 c, T2 const u) const {
    Real rho = 0., feps = 0., fsig = 0.;
    for (int n = 0; n <= NVAPOR; ++n) {
      rho += u[n];
      c[n] = u[n] / mu_ratios_[n];
      feps += c[n];
      fsig += u[n] * cv_ratios_[n];
    }
    Real KE = 0.5 * (u[IM1] * u[IM1] + u[IM2] * u[IM2] + u[IM3] * u[IM3]) / rho;
    Real gm1 = pmy_block->peos->GetGamma() - 1.;
    c[IPR] = gm1 * (u[IEN] - KE) * feps / fsig;
    c[IDN] = c[IPR] / (feps * Rd_);
    c[IVX] = u[IVX] / rho;
    c[IVY] = u[IVY] / rho;
    c[IVZ] = u[IVZ] / rho;

    Real mols = c[IPR] / (Rgas * c[IDN]);
#pragma omp simd
    for (int n = 1; n <= NVAPOR; ++n) c[n] *= mols / feps;
  }

  //! Change molar mixing ratio to density
  template <typename T1, typename T2>
  void ChemicalToConserved(T1 u, T2 const c) const {
    Real sum = 1., mols = c[IPR] / (Rgas * c[IDN]);
    for (int n = 1; n <= NVAPOR; ++n) sum += c[n] / mols * (mu_ratios_[n] - 1.);
    Real rho = c[IPR] * sum / (Rd_ * c[IDN]);
    Real cvd = Rd_ / (pmy_block->peos->GetGamma() - 1.);
    u[IDN] = rho;
    u[IEN] = 0.5 * rho * (c[IVX] * c[IVX] + c[IVY] * c[IVY] + c[IVZ] * c[IVZ]);
    for (int n = 1; n <= NVAPOR; ++n) {
      u[n] = rho * c[n] / mols * mu_ratios_[n] / sum;
      u[IDN] -= u[n];
      u[IEN] += u[n] * cv_ratios_[n] * cvd * c[IDN];
    }
    u[IEN] += u[IDN] * cvd * c[IDN];
    u[IVX] = c[IVX] * rho;
    u[IVY] = c[IVY] * rho;
    u[IVZ] = c[IVZ] * rho;
  }

  // Get thermodynamic properties
  //! polytropic index $\gamma=c_p/c_v$
  template <typename T>
  Real GetGamma(T w) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += w[n] * (cv_ratios_[n] - 1.);
      feps += w[n] * (1. / mu_ratios_[n] - 1.);
    }
    return 1. + (gamma - 1.) * feps / fsig;
  }

  template <typename T>
  Real RovRd(T w) const {
    Real feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) feps += w[n] * (1. / mu_ratios_[n] - 1.);
    return feps;
  }

  //! Temperature
  template <typename T>
  Real GetTemp(T w) const {
    return w[IPR] / (w[IDN] * Rd_ * RovRd(w));
  }

  template <typename T>
  Real GetPres(T u) const {
    Real gm1 = pmy_block->peos->GetGamma() - 1;
    Real rho = 0., fsig = 0., feps = 0.;
    for (int n = 0; n <= NVAPOR; ++n) {
      rho += u[n];
      fsig += u[n] * cv_ratios_[n];
      feps += u[n] / mu_ratios_[n];
    }
    Real KE = 0.5 * (u[IM1] * u[IM1] + u[IM2] * u[IM2] + u[IM3] * u[IM3]) / rho;
    return gm1 * (u[IEN] - KE) * feps / fsig;
  }

  template <typename T>
  Real GetChi(T w) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      qsig += w[n] * (cp_ratios_[n] - 1.);
      feps += w[n] * (1. / mu_ratios_[n] - 1.);
    }
    return (gamma - 1.) / gamma * feps / qsig;
  }

  //! c_p = c_{pd}*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  template <typename T>
  Real getSpecificCp(T w) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1.;
    for (int n = 1; n <= NVAPOR; ++n) qsig += w[n] * (cp_ratios_[n] - 1.);
    return gamma / (gamma - 1.) * Rd_ * qsig;
  }

  //! c_v = c_{vd}*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
  template <typename T>
  Real getSpecificCv(T w) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1.;
    for (int n = 1; n <= NVAPOR; ++n) qsig += w[n] * (cv_ratios_[n] - 1.);
    return 1. / (gamma - 1.) * Rd_ * qsig;
  }

  //! h = c_{pd}*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  template <typename T>
  Real getSpecificEnthalpy(T w) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1.;
    for (int n = 1; n <= NVAPOR; ++n) qsig += w[n] * (cp_ratios_[n] - 1.);
    return gamma / (gamma - 1.) * Rd_ * qsig * tem[0];
  }

  template <typename T>
  Real GetMeanMolecularWeight(T w) const {
    Real feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) feps += w[n] * (mu_ratios_[n] - 1.);
    return Rgas / Rd_ * feps;
  }

  /*! Saturation surplus for vapors can be both positive and negative
   * positive value represents supersaturation \n
   * negative value represents saturation deficit
   */
  template <typename T>
  void SaturationSurplus(Real dv[], T v, VariableType vtype) const {
    Real q1[NHYDRO];
    // mass to molar mixing ratio
    if (vtype == VariableType::prim) {
      PrimitiveToChemical(q1, v);
    } else if (vtype == VariableType::cons) {
      ConservedToChemical(q1, v);
    } else {  // VariableType::chem
      for (int n = 0; n < NHYDRO; ++n) q1[n] = v[n];
    }
    // change molar density to molar mixing ratio
    Real mols = q1[IPR] / (q1[IDN] * Rgas);
    for (int n = 1; n <= NVAPOR; ++n) q1[n] /= mols;

    for (int iv = 1; iv <= NVAPOR; ++iv) {
      int nc = q1[IDN] > t3_[iv] ? iv + NVAPOR : iv + 2 * NVAPOR;
      int ic = NHYDRO - NVAPOR + nc - 1;
      Real rate = VaporCloudEquilibrium(q1, iv, ic, t3_[iv], p3_[iv], 0.,
                                        beta_[nc], delta_[nc], true);
      dv[iv] = rate / q1[iv] * v[iv];
    }
  }

 private:
  Real ftol_;
  int max_iter_;

  //! scratch array for storing variables
  // Real w1_[NHYDRO+2*NVAPOR];

  //! scratch array for storing variables
  // Real dw_[1+NVAPOR];

  // read from inputs
  //! ideal gas constant of dry air in J/kg
  Real Rd_;
  //! ratio of mean molecular weights
  Real mu_ratios_[1 + 3 * NVAPOR];
  //! ratio of specific heat capacities at constant pressure
  Real cp_ratios_[1 + 3 * NVAPOR];
  /*! dimensionless latent heat
   *$\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r} $
   */
  Real beta_[1 + 3 * NVAPOR];
  //! triple point temperature [K]
  Real t3_[1 + NVAPOR];
  //! triple point pressure [pa]
  Real p3_[1 + NVAPOR];

  // calculated quantities
  /*! latent heat in J/kg
   *$L_i^r+\Delta c_{ij}T^r == \beta_{ij}\frac{R_d}{\epsilon_i}T^r$
   */
  Real latent_[1 + 3 * NVAPOR];

  /*! dimensionless differences in specific heat capacity
   *$(c_{ij} - c_{p,i})/R_i$
   */
  Real delta_[1 + 3 * NVAPOR];

  //! ratio of specific heat capacities at constant volume
  Real cv_ratios_[1 + 3 * NVAPOR];
};

#endif
