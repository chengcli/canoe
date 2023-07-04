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

  Real mols = c[IPR] / (c[IDN] * Constants::Rgas);
#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    c[n] *= mols / sum;
  }
}

//! Change molar mixing ratio to mass mixing ratio
template <typename T1, typename T2>
void ChemicalToPrimitive(T1 w, T2 const c) const {
  // set mass mixing ratio
  Real sum = 1., mols = c[IPR] / (c[IDN] * Constants::Rgas);
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
  Real gm1 = pmy_block_->peos->GetGamma() - 1.;
  c[IPR] = gm1 * (u[IEN] - KE) * feps / fsig;
  c[IDN] = c[IPR] / (feps * Rd_);
  c[IVX] = u[IVX] / rho;
  c[IVY] = u[IVY] / rho;
  c[IVZ] = u[IVZ] / rho;

  Real mols = c[IPR] / (Constants::Rgas * c[IDN]);
#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) c[n] *= mols / feps;
}

//! Change molar mixing ratio to density
template <typename T1, typename T2>
void ChemicalToConserved(T1 u, T2 const c) const {
  Real sum = 1., mols = c[IPR] / (Constants::Rgas * c[IDN]);
  for (int n = 1; n <= NVAPOR; ++n) sum += c[n] / mols * (mu_ratios_[n] - 1.);
  Real rho = c[IPR] * sum / (Rd_ * c[IDN]);
  Real cvd = Rd_ / (pmy_block_->peos->GetGamma() - 1.);
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
