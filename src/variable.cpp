// application
#include <application/exceptions.hpp>

// canoe
#include "variable.hpp"

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

std::ostream& operator<<(std::ostream& os, Variable::Type const& type) {
  if (type == Variable::Type::MassFrac) {
    os << "Mass Fraction";
  } else if (type == Variable::Type::MassConc) {
    os << "Mass Concentration";
  } else if (type == Variable::Type::MoleFrac) {
    os << "Mole Fraction";
  } else if (type == Variable::Type::MoleConc) {
    os << "Mole Concentration";
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, Variable const& var) {
  os << var.mytype_ << ": ";
  for (auto& v : var.data_) os << v << ", ";
  return os;
}

void Variable::ConvertTo(Variable::Type type) {
  if (type == mytype_) {
    return;
  }

  if (type == Type::MassFrac) {
    ConvertToMassFraction();
  } else if (type == Type::MassConc) {
    ConvertToMassConcentration();
  } else if (type == Type::MoleFrac) {
    ConvertToMoleFraction();
  } else if (type == Type::MoleConc) {
    ConvertToMoleConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

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
  Real g = 1.;

  // set molar mixing ratio
  Real sum = 1.;
#pragma omp simd reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }

#pragma omp simd reduction(+ : sum, g)
  for (int n = 0; n < NCLOUD; ++n) {
    sum += c[n] * (pthermo->GetInvMuRatio(n + 1 + NVAPOR) - 1.);
    g += -c[n];
  }

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetInvMuRatio(n) / sum;
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    c[n] *= pthermo->GetInvMuRatio(n + 1 + NVAPOR) / sum;
  }

  // set temperature
  w[IDN] = w[IPR] / (w[IDN] * g * pthermo->GetRd() * sum);

  SetType(Type::MoleFrac);
}

void Variable::moleFractionToMassFraction() {
  auto pthermo = Thermodynamics::GetInstance();
  Real g = 1.;

  // set mass mixing ratio
  Real sum = 1.;
#pragma omp simd reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum += w[n] * (pthermo->GetMuRatio(n) - 1.);
  }

#pragma omp simd reduction(+ : sum, g)
  for (int n = 0; n < NCLOUD; ++n) {
    sum += c[n] * (pthermo->GetMuRatio(n + 1 + NVAPOR) - 1.);
    g += -c[n];
  }

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetMuRatio(n) / sum;
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    c[n] *= pthermo->GetMuRatio(n + 1 + NVAPOR) / sum;
  }

  // set density
  w[IDN] = sum * w[IPR] / (w[IDN] * g * pthermo->GetRd());

  SetType(Type::MassFrac);
}

void Variable::massConcentrationToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  Real rho = 0., sum = 0., cvt = 0., rhoR = 0.;

#pragma omp simd reduction(+ : rho, sum, cvt, rhoR)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += w[n];
    sum += w[n] * pthermo->GetInvMuRatio(n);
    cvt += w[n] * pthermo->GetCvMassRef(n);
    rhoR += w[n] * Constants::Rgas * pthermo->GetInvMu(n);
  }

#pragma omp simd reduction(+ : rho, sum, cvt)
  for (int n = 0; n < NCLOUD; ++n) {
    rho += c[n];
    sum += c[n] * pthermo->GetInvMuRatio(n + 1 + NVAPOR);
    cvt += c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR);
  }

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= pthermo->GetInvMuRatio(n) / sum;
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    c[n] *= pthermo->GetInvMuRatio(n + 1 + NVAPOR) / sum;
  }

  Real di = 1. / rho;

  w[IDN] = w[IEN] / cvt;
  w[IPR] = rhoR * w[IDN];
  w[IVX] *= di;
  w[IVY] *= di;
  w[IVZ] *= di;

  SetType(Type::MoleFrac);
}

void Variable::moleFractionToMassConcentration() {
  auto pthermo = Thermodynamics::GetInstance();

  Real tem = w[IDN];
  Real sum = 1., g = 1.;

#pragma omp simd reduction(+ : sum)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum += w[n] * (pthermo->GetMuRatio(n) - 1.);
  }

#pragma omp simd reduction(+ : sum, g)
  for (int n = 0; n < NCLOUD; ++n) {
    sum += c[n] * (pthermo->GetMuRatio(n + 1 + NVAPOR) - 1.);
    g += -c[n];
  }

  Real rho = w[IPR] * sum / (pthermo->GetRd() * tem * g);

  w[IDN] = rho;
  w[IEN] = 0.;

  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= rho * pthermo->GetMuRatio(n) / sum;
    w[IDN] -= w[n];
    w[IEN] += w[n] * pthermo->GetCvMassRef(n) * tem;
  }

  for (int n = 0; n < NCLOUD; ++n) {
    c[n] *= rho * pthermo->GetMuRatio(n + 1 + NVAPOR) / sum;
    w[IDN] -= c[n];
    w[IEN] += c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR) * tem;
  }

  w[IEN] += w[IDN] * pthermo->GetCvMassRef(0) * tem;
  w[IVX] *= rho;
  w[IVY] *= rho;
  w[IVZ] *= rho;

  SetType(Type::MassConc);
}

void Variable::massFractionToMassConcentration() {
  auto pthermo = Thermodynamics::GetInstance();
  Real igm1 = 1.0 / (pthermo->GetGammadRef() - 1.0);

  // density
  Real rho = w[IDN], pres = w[IPR];
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] *= rho;
    w[IDN] -= w[n];
  }

  for (int n = 0; n < NCLOUD; ++n) {
    c[n] *= rho;
    w[IDN] -= c[n];
  }

  Real fsig = 1., feps = 1.;
  // vapors
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * (pthermo->GetCvRatioMass(n) - 1.);
    feps += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }

  // clouds
#pragma omp simd reduction(+ : fsig)
  for (int n = 0; n < NCLOUD; ++n) {
    fsig += c[n] * (pthermo->GetCvRatioMass(n + 1 + NVAPOR) - 1.);
  }

  // internal energy
  w[IEN] = igm1 * pres * fsig / feps;

  // momentum
  w[IVX] *= rho;
  w[IVY] *= rho;
  w[IVZ] *= rho;

  SetType(Type::MassConc);
}

void Variable::massConcentrationToMassFraction() {
  auto pthermo = Thermodynamics::GetInstance();
  Real gm1 = pthermo->GetGammadRef() - 1.;

  Real rho = 0.;
#pragma omp simd reduction(+ : rho)
  for (int n = 0; n <= NVAPOR; ++n) rho += w[n];

#pragma omp simd reduction(+ : rho)
  for (int n = 0; n < NCLOUD; ++n) rho += c[n];

  w[IDN] = rho;
  Real di = 1. / rho;

  // mass mixing ratio
#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) w[n] *= di;

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) c[n] *= di;

  Real fsig = 1., feps = 1.;
  // vapors
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * (pthermo->GetCvRatioMass(n) - 1.);
    feps += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }

#pragma omp simd reduction(+ : fsig)
  for (int n = 0; n < NCLOUD; ++n) {
    fsig += c[n] * (pthermo->GetCvRatioMass(n + 1 + NVAPOR) - 1.);
  }

  w[IPR] = gm1 * w[IEN] * feps / fsig;

  // velocity
  w[IVX] *= di;
  w[IVY] *= di;
  w[IVZ] *= di;

  SetType(Type::MassFrac);
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
