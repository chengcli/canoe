#pragma once

#include "thermodynamics.hpp"

template <typename T>
T potential_temperature(Thermodynamics const *pthermo, T w, float p0) {
  return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
}

template <typename T>
T equivalent_potential_temperature(Thermodynamics *pthermo, T w, float p0) {
  return GetTemp(w) * pow(p0 / w[IPR], GetChi(w));
}

template <typename A, typename B>
A moist_static_energy(Thermodynamics *pthermo, A w, B gz) {
  return pthermo->MoistEnthalpy(w) + gz;
}

template <typename T>
T relative_humidity(Thermodynamics *pthermo, T w, int ivapor) {
  return pthermo->GetDensity(iH2O) / w[IDN];
}

/*
template <typename T>
Real Thermodynamics::EquivalentPotentialTemp(T w) {
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
  return PotentialTemp(w, p0);
#endif
}*/
