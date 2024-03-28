// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

Real __attribute__((weak))
get_gammad(Thermodynamics const* pthermo, AirParcel const& qfrac) {
  return pthermo->GetGammadRef();
}

// Eq.94 in Li2019
Real get_rovrd(Thermodynamics const* pthermo, AirParcel const& qfrac) {
  qfrac.ToMoleFraction();
  Real fgas = 1., feps = 1.;

#pragma omp simd reduction(+ : feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    feps += qfrac.w[n] * (pthermo->GetMuRatio(n) - 1.);
  }

#pragma omp simd reduction(+ : fgas, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    fgas += -qfrac.c[n];
    feps += qfrac.c[n] * (pthermo->GetMuRatio(n) - 1.);
  }

  return fgas / feps;
}

// TODO(cli): check
Real get_chi(Thermodynamics const* pthermo, AirParcel const& qfrac) {
  qfrac.ToMoleFraction();
  Real gammad = get_gammad(pthermo, qfrac);

  Real qsig = 1., feps = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n) {
    qsig += qfrac.w[n] * (pthermo->GetCpRatioMole(n) - 1.);
  }

#pragma omp simd reduction(+ : qsig, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    feps += -qfrac.c[n];
    qsig += qfrac.c[n] * (pthermo->GetCpRatioMole(n + 1 + NVAPOR) - 1.);
  }

  return (gammad - 1.) / gammad / qsig;
}

//! Eq.XX in Li2019
Real cal_dlnT_dlnP(Thermodynamics const* pthermo, AirParcel const& qfrac,
                   Real latent[]) {
  qfrac.ToMoleFraction();
  // calculate gammad
  Real gammad = get_gammad(pthermo, qfrac);

  Real q_gas = 1.;
  for (int n = 0; n < NCLOUD; ++n) q_gas -= qfrac.c[n];

  Real f_sig = 1.;
  // vapor
  for (int n = 1; n <= NVAPOR; ++n)
    f_sig += qfrac.w[n] * (pthermo->GetCpRatioMole(n) - 1.);
  // cloud
  for (int n = 0; n < NCLOUD; ++n)
    f_sig += qfrac.c[n] * (pthermo->GetCpRatioMole(1 + NVAPOR + n) - 1.);
  Real cphat_ov_r = gammad / (gammad - 1.) * f_sig / q_gas;

  // vapor
  Real xd = q_gas;
  for (int n = 1; n <= NVAPOR; ++n) xd -= qfrac.w[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    c1 += qfrac.w[iv] / xd * latent[iv];
    c2 += qfrac.w[iv] / xd * latent[iv] * latent[iv];
    c3 += qfrac.w[iv] / xd;
  }

  return (1. + c1) / (cphat_ov_r + (c2 + c1 * c1) / (1. + c3));
}

void set_total_equivalent_vapor(
    AirParcel* qfrac, std::vector<IndexSet> const& cloud_index_set,
    std::map<IndexPair, ReactionInfo> const& cloud_reaction_map) {
  // vpaor <=> cloud
  for (int i = 1; i <= NVAPOR; ++i)
    for (auto& j : cloud_index_set[i]) {
      qfrac->w[i] += qfrac->c[j];
      qfrac->c[j] = 0.;
    }

  // vapor + vapor <=> cloud
  for (auto const& [ij, info] : cloud_reaction_map) {
    auto& indx = info.first;
    auto& stoi = info.second;

    qfrac->w[indx[0]] += stoi[0] / stoi[2] * qfrac->c[indx[2]];
    qfrac->c[indx[2]] = 0.;
    qfrac->w[indx[1]] += stoi[1] / stoi[2] * qfrac->c[indx[2]];
    qfrac->c[indx[2]] = 0.;
  }
}
