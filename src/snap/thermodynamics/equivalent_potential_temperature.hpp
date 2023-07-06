
//! Equivalent potential temperature
/*Real MoistEntropy(AthenaArray<Real> const& w, Thermodynamics *pthermo,
Particles *ppart, int k, int j, int i) { #if (NVAPOR > 0) Real gamma =
pthermo->pmy_block->peos->GetGamma(); Real tem[1] = {pthermo->GetTemp(prim)};
  update_gamma(gamma, tem);
  Real cpd = Rd_*gamma/(gamma - 1.);
  Real temp = tem[0]
  Real pres = prim[IPR];

  Real qd = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qd -= prim[n];

  Real qc[1+NVAPOR];
  std::fill(qc, qc + 1 + NVAPOR, 0.);
  for (int n = 1 + NVAPOR; n <= NVAPOR; ++n)
    qc[1+(n-1)%NVAPOR] += prim[n] + 1.0E-10;  // prevent devide by 0

  Real lv = 0.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    lv += GetLatent(n,temp)*prim[ng]*ratio;
  }

  Real st = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    st += (prim[n] + prim[ng]*ratio)*(GetCpRatio(n) - 1.);
  }
  Real lv_ov_cpt = lv/(cpd*st*temp);

  Real chi = Rd_/cpd*qd/st;

  Real xv = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    xv += prim[n]/qd/mu_ratios_[n];

  Real pd = pres/xv;

  Real rh = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    Real eta = prim[n]/qd/mu_ratios_[n];
    Real pv = pres*eta/xv;
    int nc = n + NVAPOR;
    Real esat;
    if (n == AMMONIA_VAPOR_ID)
      esat = sat_vapor_p_NH3_BriggsS(temp);
    else if (n == WATER_VAPOR_ID)
      esat = sat_vapor_p_H2O_BriggsS(temp);
    else
      esat = SatVaporPresIdeal(temp/t3_[nc], p3_[nc], beta_[nc], delta_[nc]);
    rh *= pow(pv/esat, -eta*Rd_/(cpd*st));
  }

  return temp*pow(p0/pd, chi)*exp(lv_ov_cpt)*rh;
#else
  return PotentialTemp(prim, p0, pthermo);
#endif
}*/
