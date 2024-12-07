// C/C++
#include <algorithm>
#include <iomanip>
#include <iostream>

// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::_rk4_integrate_z(Real dz, std::string method, Real grav,
                                      Real adTdz) const {
  auto& thermo = kinetics_->thermo();

  Real step[] = {0.5, 0.5, 1.};
  Real dTdz[4], chi[4];
  Real latent[1 + NVAPOR];
  Real cp_ratio_mole[Size];

  std::array<Real, Size> enthalpy;
  std::array<Real, Size> xfrac;

  Real temp0 = thermo.temperature();
  Real pres0 = thermo.pressure();

  for (int i = 0; i < Size; ++i) {
    cp_ratio_mole[i] = cp_ratio_[i] / inv_mu_ratio_[i];
  }

  for (int rk = 0; rk < 4; ++rk) {
    thermo.getMoleFractions(xfrac.data());
    thermo.getEnthalpy_RT(enthalpy.data());

    for (int i = 1; i <= NVAPOR; ++i) {
      if (xfrac[i + NVAPOR] > 0.) {
        latent[i] = enthalpy[i] - enthalpy[i + NVAPOR];
      } else {
        latent[i] = 0.;
      }
    }

    if (method != "reversible") {
      for (int j = 1 + NVAPOR; j < Size; ++j) xfrac[j] = 0;
    }

    Real q_gas = 1., q_eps = 1.;
    for (int i = 1; i <= NVAPOR; ++i) {
      q_eps += xfrac[i] * (1. / inv_mu_ratio_[i] - 1.);
    }

    for (int j = 1 + NVAPOR; j < Size; ++j) {
      q_eps += xfrac[j] * (1. / inv_mu_ratio_[j] - 1.);
      q_gas += -xfrac[j];
    }

    Real g_ov_Rd = grav / Rd_;
    Real R_ov_Rd = q_gas / q_eps;

    if (method == "reversible" || method == "pseudo") {
      chi[rk] = cal_dlnT_dlnP(xfrac.data(), gammad_, cp_ratio_mole, latent);
    } else if (method == "dry") {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = cal_dlnT_dlnP(xfrac.data(), gammad_, cp_ratio_mole, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }

    dTdz[rk] = -chi[rk] * g_ov_Rd / R_ov_Rd + adTdz;
    chi[rk] = -R_ov_Rd / g_ov_Rd * dTdz[rk];

    // integrate over dz
    Real chi_avg, temp, pres;
    if (rk < 3) {
      temp = temp0 + dTdz[rk] * dz * step[rk];
      chi_avg = chi[rk];
    } else {
      temp = temp0 +
             1. / 6. * (dTdz[0] + 2. * dTdz[1] + 2. * dTdz[2] + dTdz[3]) * dz;
      chi_avg = 1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]);
    }

    if (!(temp > 0.)) temp = temp0;

    if (fabs(temp - temp0) > 0.01) {
      pres = pres0 * pow(temp / temp0, 1. / chi_avg);
    } else {  // isothermal limit
      pres = pres0 * exp(-2. * g_ov_Rd * dz / (R_ov_Rd * (temp + temp0)));
    }

    thermo.setMoleFractions(xfrac.data());
    EquilibrateTP(temp, pres);
  }
}
