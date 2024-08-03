// C/C++
#include <algorithm>
#include <iostream>

// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::_rk4_integrate_lnp(Real dlnp, std::string method,
                                        Real adlnTdlnP) const {
  auto& thermo = kinetics_->thermo();

  Real step[] = {0.5, 0.5, 1.};
  Real chi[4];
  Real latent[1 + NVAPOR];
  Real cp_ratio_mole[Size];

  std::vector<Real> enthalpy(Size);
  std::vector<Real> xfrac(Size);

  Real temp0 = thermo.temperature();
  Real pres0 = thermo.pressure();

  for (int i = 0; i < Size; ++i) {
    cp_ratio_mole[i] = cp_ratio_[i] / inv_mu_ratio_[i];
  }

  for (int rk = 0; rk < 4; ++rk) {
    thermo.getEnthalpy_RT(enthalpy.data());
    thermo.getMoleFractions(xfrac.data());

    for (int i = 1; i <= NVAPOR; ++i) {
      if (thermo.massFraction(i + NVAPOR) > 0.) {
        latent[i] = enthalpy[i] - enthalpy[i + NVAPOR];
      } else {
        latent[i] = 0.;
      }
      latent[i] *= temp0 * Constants::Rgas;
    }

    if (method != "reversible") {
      for (int j = 1 + NVAPOR; j < Size; ++j) xfrac[j] = 0;
    }

    // calculate tendency
    if (method == "reversible" || method == "pseudo") {
      chi[rk] = cal_dlnT_dlnP(xfrac.data(), gammad_, cp_ratio_mole, latent);
    } else if (method == "dry") {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = cal_dlnT_dlnP(xfrac.data(), gammad_, cp_ratio_mole, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }
    chi[rk] += adlnTdlnP;

    Real temp, pres;
    // integrate over dlnp
    if (rk < 3) {
      temp = temp0 * exp(chi[rk] * dlnp * step[rk]);
    } else {
      temp = temp0 * exp(1. / 6. *
                         (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp);
    }
    pres = pres0 * exp(dlnp);

    thermo.setMoleFractions(xfrac.data());
    EquilibrateTP(temp, pres);
  }
}
