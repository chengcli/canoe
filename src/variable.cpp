// application
#include <application/exceptions.hpp>

// canoe
#include "variable.hpp"

// climath
#include <climath/core.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

void Variable::ConvertToMassFraction() {
  if (mytype_ == Type::MassFrac) {
    return;
  }

  if (mytype_ == Type::MassConc) {
    massConcentrationToMassFraction();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToMassFraction();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToMassFraction();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

void Variable::ConvertToMassConcentration() {
  if (mytype_ == Type::MassConc) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    massFractionToMassConcentration();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToMassConcentration();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToMassConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

void Variable::ConvertToMoleFraction() {
  if (mytype_ == Type::MoleFrac) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    massFractionToMoleFraction();
  } else if (mytype_ == Type::MassConc) {
    massConcentrationToMoleFraction();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToMoleFraction();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

void Variable::ConvertToMoleConcentration() {
  if (mytype_ == Type::MoleConc) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    massFractionToMoleConcentration();
  } else if (mytype_ == Type::MassConc) {
    massConcentrationToMoleConcentration();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToMoleConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

void Variable::massFractionToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  // set molar mixing ratio
  Real sum = 1.;
#pragma reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetInvMuRatio(n) / sum;
  }

  // set temperature
  w[IDN] = w[IPR] / (w[IDN] * pthermo->GetRd() * sum);

  mytype_ = Type::MoleFrac;
}

void Variable::moleFractionToMassFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  // set mass mixing ratio
  Real sum = 1.;
#pragma reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum += w[n] * (pthermo->GetMuRatio(n) - 1.);
  }

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetMuRatio(n) / sum;
  }

  // set density
  w[IDN] = sum * w[IPR] / (w[IDN] * pthermo->GetRd());

  mytype_ = Type::MassFrac;
}

void Variable::massConcentrationToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  Real rho = 0., feps = 0., fsig = 0.;

#pragma reduction(+ : rho, feps, fsig)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += w[n];
    feps += w[n] * pthermo->GetInvMuRatio(n);
    fsig += w[n] * pthermo->GetCvRatioMass(n);
  }

#pragma omp simd
  for (int n = 0; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetInvMuRatio(n);
  }

  Real di = 1. / rho;

  // TODO(cli): not true for cubed-sphere
  Real KE = 0.5 * (sqr(w[IM1]) + sqr(w[IM2]) + sqr(w[IM3])) * di;
  Real gm1 = pthermo->GetGammad(*this) - 1.;

  w[IPR] = gm1 * (w[IEN] - KE) * feps / fsig;
  w[IDN] = w[IPR] / (feps * pthermo->GetRd());
  w[IVX] *= di;
  w[IVY] *= di;
  w[IVZ] *= di;

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) w[n] /= feps;

  mytype_ = Type::MassFrac;
}

void Variable::moleFractionToMassConcentration() {
  auto pthermo = Thermodynamics::GetInstance();

  Real tem = w[IDN];
  Real sum = 1.;

#pragma reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) sum += w[n] * (pthermo->GetMuRatio(n) - 1.);

  Real rho = w[IPR] * sum / (pthermo->GetRd() * tem);

  w[IDN] = rho;

  // TODO(cli): not true for cubed-sphere
  w[IEN] = 0.5 * rho * (sqr(w[IVX]) + sqr(w[IVY]) + sqr(w[IVZ]));

  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= rho * pthermo->GetMuRatio(n) / sum;
    w[IDN] -= w[n];
    w[IEN] += w[n] * pthermo->GetCvMass(*this, n) * tem;
  }

  w[IEN] += w[IDN] * pthermo->GetCvMass(*this, 0) * tem;
  w[IVX] *= rho;
  w[IVY] *= rho;
  w[IVZ] *= rho;

  mytype_ = Type::MassConc;
}

void Variable::massFractionToMassConcentration() {
  auto pthermo = Thermodynamics::GetInstance();
  Real igm1 = 1.0 / (pthermo->GetGammad(*this) - 1.0);

  // density
  Real rho = w[IDN], pres = w[IPR];
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= rho;
    w[IDN] -= w[n];
  }

  // total energy
  w[IEN] = 0.5 * rho * (sqr(w[IVX]) + sqr(w[IVY]) + sqr(w[IVZ]));
  Real fsig = 1., feps = 1.;
  // vapors
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * (pthermo->GetCvRatioMass(n) - 1.);
    feps += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }
  w[IEN] += igm1 * pres * fsig / feps;

  // momentum
  w[IVX] *= rho;
  w[IVY] *= rho;
  w[IVZ] *= rho;

  mytype_ = Type::MassConc;
}

void Variable::massConcentrationToMassFraction() {
  auto pthermo = Thermodynamics::GetInstance();
  Real gm1 = pthermo->GetGammad(*this) - 1.;

  Real rho = 0.;
  for (int n = 0; n <= NVAPOR; ++n) rho += w[n];
  w[IDN] = rho;
  Real di = 1. / rho;

  // mass mixing ratio
#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) w[n] *= di;

  // pressure
  Real KE = 0.5 * di * (sqr(w[IVX]) + sqr(w[IVY]) + sqr(w[IVZ]));
  Real fsig = 1., feps = 1.;
  // vapors
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * (pthermo->GetCvRatioMass(n) - 1.);
    feps += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }
  w[IPR] = gm1 * (w[IEN] - KE) * feps / fsig;

  // velocity
  w[IVX] *= di;
  w[IVY] *= di;
  w[IVZ] *= di;

  mytype_ = Type::MassFrac;
}

void Variable::moleFractionToMoleConcentration() {
  throw NotImplementedError("Variable::moleFractionToMoleConcentration");
}

void Variable::moleConcentrationToMoleFraction() {
  throw NotImplementedError("Variable::moleConcentrationToMoleFraction");
}

void Variable::massFractionToMoleConcentration() {
  throw NotImplementedError("Variable::massFractionToMoleConcentration");
}

void Variable::moleConcentrationToMassFraction() {
  throw NotImplementedError("Variable::moleConcentrationToMassFraction");
}

void Variable::massConcentrationToMoleConcentration() {
  throw NotImplementedError("Variable::massConcentrationToMoleConcentration");
}

void Variable::moleConcentrationToMassConcentration() {
  throw NotImplementedError("Variable::moleConcentrationToMassConcentration");
}
