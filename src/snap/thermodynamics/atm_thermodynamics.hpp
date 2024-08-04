#pragma once

#include "thermodynamics.hpp"

template <typename T>
Real potential_temp(Thermodynamics const *pthermo, T w, Real p0) {
  return pthermo->GetTemp(w) * pow(p0 / w[IPR], pthermo->GetChi(w));
}

template <typename T>
Real moist_static_energy(Thermodynamics const *pthermo, T w, Real gz) {
  pthermo->SetPrimitive(w);
  Real intEng = pthermo->GetInternalEnergy(w);
  Real tempv = pthermo->GetTemp() * pthermo->RovRd();
  return intEng + pthermo->GetRd() * tempv + gz;
}

template <typename T>
std::array<Real, Thermodynamics::Size> relative_humidity(
    Thermodynamics const *pthermo, T w) {
  pthermo->SetPrimitive(w);
  auto kinetics = get_kinetics_object(pthermo);

  std::vector<Real> kfwd(kinetics->nReactions());
  kinetics->getFwdRateConstants(kfwd.data());

  std::array<Real, Thermodynamics::Size> xfrac;
  kinetics->thermo().getMoleFractions(xfrac.data());

  // gas fractions
  Real xg = 0.;
  for (int n = 0; n <= NVAPOR; ++n) xg += xfrac[n];

  Real temp = pthermo->GetTemp();
  for (int n = 1; n <= NVAPOR; ++n) {
    xfrac[n] *= w[IPR] / (xg * kfwd[n - 1] * Cantera::GasConstant * temp);
  }

  return xfrac;
}

// Eq.4.5.11 in Emanuel (1994)
template <typename T>
Real equivalent_potential_temp(Thermodynamics const *pthermo, T w, Real rh,
                               Real p0) {
  pthermo->SetPrimitive(w);
  auto kinetics = get_kinetics_object(pthermo);
  auto &thermo = kinetics->thermo();

  std::array<Real, Thermodynamics::Size> enthalpy;
  std::array<Real, Thermodynamics::Size> cp;
  thermo.getEnthalpy_RT(enthalpy.data());
  thermo.getPartialMolarCp(cp.data());

  Real temp = thermo.temperature();
  Real pres = thermo.pressure();

  Real lv = (enthalpy[1] - enthalpy[2]) * Cantera::GasConstant * temp /
            thermo.molecularWeight(1);
  Real cpd = cp[0] / thermo.molecularWeight(0);
  Real cl = cp[2] / thermo.molecularWeight(2);

  /*std::cout << "lv: " << lv << std::endl;
  std::cout << "cpd: " << cpd << std::endl;
  std::cout << "cl: " << cl << std::endl;*/

  Real qd = thermo.massFraction(0);
  Real qv = thermo.massFraction(1);
  Real qc = thermo.massFraction(2);
  Real qt = qv + qc;

  Real xd = thermo.moleFraction(0);
  Real xv = thermo.moleFraction(1);
  Real xg = xd + xv;

  Real Rd = Cantera::GasConstant / thermo.molecularWeight(0);
  Real Rv = Cantera::GasConstant / thermo.molecularWeight(1);
  Real pd = xd / xg * pres;
  Real cpt = cpd * qd + cl * qt;

  /*std::cout << "Rd: " << Rd << std::endl;
  std::cout << "Rv: " << Rv << std::endl;
  std::cout << "pd: " << pd << std::endl;
  std::cout << "cpt: " << cpt << std::endl;
  std::cout << "rh: " << rh << std::endl;*/

  return temp * pow(p0 / pd, Rd * qd / cpt) * pow(rh, -Rv * qv / cpt) *
         exp(lv * qv / (cpt * temp));
}
