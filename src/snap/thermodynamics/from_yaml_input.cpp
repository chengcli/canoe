
void EquationOfState::InitFromYAMLInput(std::string const& fname) {
  auto thermo = Cantera::newThermo(fname, "gas");

  if (thermo->nSpecies() != IVX) {
    throw RuntimeError("Thermodynamics",
                       "Number of species does not match the input file");
  }

  kinetics_ = Cantera::newKinetics({thermo}, fname);
  UpdateThermoProperty();
}

void EquationOfState::UpdateThermoProperty() {
  auto thermo = kinetics_->thermo();

  // --------- vapor + cloud thermo ---------
  std::vector<Real> mu(thermo->nSpecies());
  std::vector<Real> cp_mole(thermo->nSpecies());
  std::vector<Real> cp_ratio_mole(thermo->nSpecies());

  // g/mol
  thermo->getMolecularWeights(mu.data());

  // J/kmol/K
  thermo->getPartialMolarCp(cp_mole.data());

  Rd_ = Cantera::GasConstant / mu[0];
  gammad_ = cp_mole[0] / (cp_mole[0] - Cantera::GasConstant);

  // ---------- dimensionless properties ----------

  // cp ratios
  for (size_t i = 0; i < cp_mole.size(); ++i) {
    cp_ratio_mole[i] = cp_mole[i] / cp_mole[0];
  }

  // molecular weight ratios
  for (size_t i = 0; i < mu.size(); ++i) {
    inv_mu_ratio_[i] = mu[0] / mu[i];
    cp_ratio_mass[n] = cp_ratio_mole[n] * inv_mu_ratio[n];
  }

  // set up cv ratio
  // calculate cv_ratio = $\sigma_i + (1. - \gamma)/\epsilon_i$
  for (int n = 0; n <= NVAPOR; ++n) {
    cv_ratio_mass_[n] =
        gammad * cp_ratio_mass[n] + (1. - gammad) * inv_mu_ratio[n];
  }

  for (int n = 1 + NVAPOR; n < IVX; ++n) {
    cv_ratio_mass_[n] = gammad * cp_ratio_mass[n];
  }
}
