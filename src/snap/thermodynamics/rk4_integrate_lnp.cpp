// C/C++
#include <algorithm>
#include <iostream>

// thermodynamics
#include "thermodynamics.hpp"

void Thermodynamics::_rk4_integrate_lnp(Real dlnp, std::string method,
                                        Real adlnTdlnP) {
  auto thermo = kinetics_->thermo();

  Real step[] = {0.5, 0.5, 1.};
  Real chi[4];
  Real latent[1 + NVAPOR];
  Real cp_ratio_mole[Size];

  std::vector<Real> rates(Size);
  std::vector<Real> enthalpy(Size);
  std::vector<Real> xfrac(Size);

  thermo.getEnthalpy_RT(enthalpy.data());
  kinetics_->getNetProductionRates(rates.data());

  for (int i = 1; i <= NVAPOR; ++i) {
    if (rates[i] > 0.) {
      latent[i] = (enthalpy[i] - enthalpy[i + NVAPOR]) * thermo.temperature() *
                  Constants::Rgas;
    } else {
      latent[i] = 0.;
    }
  }

  for (int i = 0; i < Size; ++i) {
    cp_ratio_mole[i] = cp_ratio_[i] / inv_mu_ratio_[i];
  }

  for (int rk = 0; rk < 4; ++rk) {
    EquilibrateTP();
    thermo.getMoleFractions(xfrac.data());

    Real temp = thermo.temperature();
    Real pres = thermo.pressure();

    if (method != "reversible") {
      for (int j = 1 + NVAPOR; j < Size; ++j) xfrac[j] = 0;
    }
    thermo.setMoleFractions(xfrac.data());
    thermo.setPressure(pres);

    // calculate tendency
    if (method == "reversible" || method == "pseudo") {
      chi[rk] = cal_dlnT_dlnP(xfrac, gammad_, cp_ratio_mole, latent);
    } else if (method == "dry") {
      for (int i = 1; i <= NVAPOR; ++i) latent[i] = 0;
      chi[rk] = cal_dlnT_dlnP(xfrac, gammad_, cp_ratio_mole, latent);
    } else {  // isothermal
      chi[rk] = 0.;
    }
    chi[rk] += adlnTdlnP;

    // integrate over dlnp
    if (rk < 3) {
      thermo.setTemperature(temp * exp(chi[rk] * dlnp * step[rk]));
    } else {
      thermo.setTemperature(
          temp *
          exp(1. / 6. * (chi[0] + 2. * chi[1] + 2. * chi[2] + chi[3]) * dlnp));
    }
    thermo.setPressure(pres * exp(dlnp));
  }

  // recondensation
  EquilibrateTP();
  thermo.getMoleFractions(xfrac.data());
  pres = thermo.pressure();
  if (method != "reversible") {
    for (int j = 1 + NVAPOR; j < Size; ++j) xfrac[j] = 0;
  }
  thermo.setMoleFractions(xfrac.data());
  thermo.setPressure(pres);
}
