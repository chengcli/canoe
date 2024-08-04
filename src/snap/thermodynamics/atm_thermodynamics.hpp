#pragma once

#include "thermodynamics.hpp"

template <typename A, typename B>
B potential_temp(Thermodynamics const *pthermo, A w, B p0) {
  return pthermo->GetTemp(w) * pow(p0 / w[IPR], pthermo->GetChi(w));
}

template <typename A, typename B>
B moist_static_energy(Thermodynamics const *pthermo, A w, B gz) {
  pthermo->SetPrimitive(w);
  Real intEng = pthermo->GetInternalEnergy(w);
  Real tempv = pthermo->GetTemp() * pthermo->RovRd();
  return intEng + pthermo->GetRd() * tempv + gz;
}

template <typename T>
std::vector<Real> relative_humidity(Thermodynamics const *pthermo, T w) {
  pthermo->SetPrimitive(w);
  auto kinetics = get_kinetics_object(pthermo);

  std::vector<Real> kfwd(kinetics->nReactions());
  kinetics->getFwdRateConstants(kfwd.data());

  std::vector<Real> xfrac(Thermodynamics::Size);
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

/*template <typename T>
Real Thermodynamics::equivalent_potential_temp(Thermodynamics const *pthermo,
    T w) {
#if (NVAPOR > 0)
  Real cpd = GetRd() * gammad_ / (gammad_ - 1.);
  Real temp = GetTemp(w);
  Real pres = w[IPR];

  Real qd = 1.;
  for (int n = 1; n < Size; ++n) qd -= w[n];

  Real qc[1 + NVAPOR];
  std::fill(qc, qc + 1 + NVAPOR, 0.);
  for (int n = 1 + NVAPOR; n < Size; ++n)
    qc[1 + (n - 1) % NVAPOR] += w[n] + 1.0E-10;  // prevent devide by 0

  Real lv = 0.;
  for (int n = 1 + NVAPOR; n < Size; ++n) {
    int ng = 1 + (n - 1) % NVAPOR;
    Real ratio = (w[n] + 1.0E-10) / qc[ng];
    lv += GetLatent_RT(n - NVAPOR, n) * w[ng] * ratio;
  }

  Real st = 1.;
  for (int n = 1 + NVAPOR; n < Size; ++n) {
    int ng = 1 + (n - 1) % NVAPOR;
    Real ratio = (w[n] + 1.0E-10) / qc[ng];
    st += (w[n] + w[ng] * ratio) * (cp_ratio_[n] - 1.);
  }
  Real lv_ov_cpt = lv / (cpd * st * temp);

  Real chi = GetRd() / cpd * qd / st;

  Real xv = 1.;
  for (int n = 1; n <= NVAPOR; ++n) xv += w[n] / qd * inv_mu_ratio_[n];

  Real pd = pres / xv;
  Real rh = 1.;

  std::array<Real, Size> kfwd;
  kinetics_->getFwdRateConstants(kfwd.data());

  for (int n = 1; n <= NVAPOR; ++n) {
    Real eta = w[n] / qd * inv_mu_ratio_[n];
    Real pv = pres * eta / xv;
    int nc = n + NVAPOR;
    Real esat = std::max(pv, kfwd[j] * Cantera::GasConstant * temp);
    rh *= pow(pv / esat, -eta * GetRd() / (cpd * st));
  }

  return temp * pow(p0 / pd, chi) * exp(lv_ov_cpt) * rh;
#else
  return potential_temp(pthermo, w, p0);
#endif
}*/
