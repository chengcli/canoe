#include "thermodynamics.hpp"

template <typename T>
Real Thermodynamics::RovRd(T w) const {
  Real feps = 1.;
  for (int n = 1 + NVAPOR; n <= Size; ++n) feps -= w[n];

  for (int n = 1; n <= NVAPOR; ++n) feps += w[n] * (inv_mu_ratio_[n] - 1.);
  return feps;
}

template <typename T>
Real Thermodynamics::GetChi(T w) const {
  Real qsig = 1., feps = 1.;
  for (int n = 1 + NVAPOR; n <= Size; ++n) {
    feps -= w[n];
    qsig += w[n] * (cp_ratio_[n] - 1.);
  }

  for (int n = 1; n <= NVAPOR; ++n) {
    feps += w[n] * (inv_mu_ratio_[n] - 1.);
    qsig += w[n] * (cp_ratio_[n] - 1.);
  }

  return (gammad_ - 1.) / gammad_ * feps / qsig;
}

template <typename T>
Real Thermodynamics::GetGamma(T w) const {
  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * (cv_ratio_[n] - 1.);
    feps += w[n] * (inv_mu_ratio_[n] - 1.);
  }

  for (int n = 1 + NVAPOR; n < Size; ++n) {
    fsig += w[n] * (cv_ratio[n] - 1.);
    feps += -w[n];
  }

  return 1. + (gammad_ - 1.) * feps / fsig;
}

template <typename T>
Real Thermodynamics::GetPres(T u, std::vector<Real> cos_theta) const {
  Real igm1 = 1. / (gammad_ - 1.);

  Real rho = 0.;
  for (int n = 0; n <= NVAPOR; ++n) rho += u[n];

  // internal energy
  Real fsig = 1., feps = 1.;
  // clouds
  for (int n = 1 + NVAPOR; n < Size; ++n) {
    fsig += u[n] * (cv_ratio_[n] - 1.);
    feps -= u[n];
  }
  // vapors
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += u[n] * (cv_ratio_[n] - 1.);
    feps += u[n] * (inv_mu_ratio_[n] - 1.);
  }

  auto w = vec_raise(u, cos_theta);

  Real KE = 0.5 * (u[IM1] * w[0] + u[IM2] * w[1] + u[IM3] * w[2]) / rho;
  return igm1 * (u[IEN] - KE) * feps / fsig;
}

template <typename T>
void Thermodynamics::SetStateFromPrimitive(T w) const {
  auto thermo = kinetics_->thermo();

  thermo.setMassFractionsPartial(w);
  thermo.setDensity(w[IDN]);
  thermo.setPressure(w[IPR]);
}

template <typename T>
void Thermodynamics::SetStateFromConserved(T u) const {
  auto thermo = kinetics_->thermo();

  thermo.setMassFractions(u);
  Real rho = 0.;
  for (int n = 0; n < Size; ++n) rho += u[n];
  thermo.setDensity(rho);
  thermo.setPressure(GetPres(u));
}

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
    st += (w[n] + w[ng] * ratio) * (cp_ratio_(n) - 1.);
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
}

template <typename T>
std::vector<Real> Thermodynamics::SaturationSurplus(T w) {
  std::array<Real, 1 + NVAPOR> dq;
  Real temp = GetTemp(w);

  // mass to molar mixing ratio
  Real sum = 1.;
  for (int n = 1; n < Size; ++n) sum -= w[n];  // right now sum = qd
  dq[0] = sum;
  for (int n = 1; n <= NVAPOR; ++n) {
    dq[n] = w[n] * inv_mu_ratio_[n];
    sum += dq[n];
  }
  Real qtol = sum;
  for (int n = NVAPOR + 1; n < Size; ++n) qtol += w[n] * inv_mu_ratio_[n];

  for (int n = 0; n <= NVAPOR; ++n) dq[n] /= sum;
  // Saturation surplus for vapors can be both positive and negative
  // positive value represents supersaturation
  // negative value represents saturation deficit

  std::array<Real, Size> kfwd;
  kinetics_->getFwdRateConstants(kfwd.data());

  for (int n = 1; n <= NVAPOR; ++n) {
    Real q = dq[n];
    Real yy = q / (1. - q);
    dq[n] = -1.0E10;

    Real svp = kfwd[j] * Cantera::GasConstant * temp;
    Real xs = svp / w[IPR];
    // default to boilding (evaporate all)
    if (xs < 1.) {
      Real ys = xs / (1. - xs);
      dq[n] = std::max(dq[n], w[n] * (1. - ys / yy) * sum / qtol);
    }
  }

  return dq;
}
