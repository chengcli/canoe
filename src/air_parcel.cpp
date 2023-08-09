// application
#include <application/exceptions.hpp>

// canoe
#include "air_parcel.hpp"

// snap
#include "snap/thermodynamics/thermodynamics.hpp"

std::ostream& operator<<(std::ostream& os, AirParcel::Type const& type) {
  if (type == AirParcel::Type::MassFrac) {
    os << "Mass Fraction";
  } else if (type == AirParcel::Type::MassConc) {
    os << "Mass Concentration";
  } else if (type == AirParcel::Type::MoleFrac) {
    os << "Mole Fraction";
  } else if (type == AirParcel::Type::MoleConc) {
    os << "Mole Concentration";
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, AirParcel const& var) {
  os << var.mytype_ << ": ";
  for (auto& v : var.data_) os << v << ", ";
  return os;
}

AirParcel& AirParcel::ConvertTo(AirParcel::Type type) {
  if (type == mytype_) {
    return *this;
  }

  if (type == Type::MassFrac) {
    ToMassFraction();
  } else if (type == Type::MassConc) {
    ToMassConcentration();
  } else if (type == Type::MoleFrac) {
    ToMoleFraction();
  } else if (type == Type::MoleConc) {
    ToMoleConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  return *this;
}

AirParcel& AirParcel::ToMassFraction() {
  if (mytype_ == Type::MassFrac) {
    return *this;
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

  return *this;
}

AirParcel& AirParcel::ToMassConcentration() {
  if (mytype_ == Type::MassConc) {
    return *this;
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

  return *this;
}

AirParcel& AirParcel::ToMoleFraction() {
  if (mytype_ == Type::MoleFrac) {
    return *this;
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

  return *this;
}

AirParcel& AirParcel::ToMoleConcentration() {
  if (mytype_ == Type::MoleConc) {
    return *this;
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

  return *this;
}

void AirParcel::massFractionToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  // set molar mixing ratio
  Real sum1 = 1.;
#pragma omp simd reduction(+ : sum1)
  for (int n = 1; n <= NVAPOR; ++n) {
    sum1 += w[n] * (pthermo->GetInvMuRatio(n) - 1.);
  }

  Real sum = sum1;
#pragma omp simd reduction(+ : sum, sum1)
  for (int n = 0; n < NCLOUD; ++n) {
    sum += c[n] * (pthermo->GetInvMuRatio(n + 1 + NVAPOR) - 1.);
    sum1 += -c[n];
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
  w[IDN] = w[IPR] / (w[IDN] * pthermo->GetRd() * sum1);

  SetType(Type::MoleFrac);
}

void AirParcel::moleFractionToMassFraction() {
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

void AirParcel::massConcentrationToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();

  Real rho = 0., sum = 0., cvt = 0., rhoR = 0.;
  Real inv_rhod = 1. / w[IDN];

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

  // tracer
#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= inv_rhod;

  SetType(Type::MoleFrac);
}

void AirParcel::moleFractionToMassConcentration() {
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

  // tracer
#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= w[IDN];

  SetType(Type::MassConc);
}

void AirParcel::massFractionToMassConcentration() {
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

  // tracer
#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= w[IDN];

  SetType(Type::MassConc);
}

void AirParcel::massConcentrationToMassFraction() {
  auto pthermo = Thermodynamics::GetInstance();
  Real gm1 = pthermo->GetGammadRef() - 1.;

  Real rho = 0., inv_rhod = w[IDN];
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

  // tracer
#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= inv_rhod;

  SetType(Type::MassFrac);
}

void AirParcel::moleFractionToMoleConcentration() {
  auto pthermo = Thermodynamics::GetInstance();
  Real tem = w[IDN];
  Real pres = w[IPR];

  // total gas moles / m^3
  Real gmols = pres / (Constants::Rgas * tem);

  Real xgas = 1., fsig = 1., LE = 0.;
#pragma omp simd reduction(+ : xgas, fsig, LE)
  for (int n = 0; n < NCLOUD; ++n) {
    xgas += -c[n];
    fsig += (pthermo->GetCvRatioMole(n + 1 + NVAPOR) - 1.) * c[n];
    LE += -c[n] * pthermo->GetLatentEnergyMole(n + 1 + NVAPOR);
  }

  Real xd = xgas;
#pragma omp simd reduction(+ : xd, fsig)
  for (int n = 1; n <= NVAPOR; ++n) {
    xd += -w[n];
    fsig += (pthermo->GetCvRatioMole(n) - 1.) * w[n];
  }

  Real tmols = gmols / xgas;

  w[IDN] = tmols * xd;
  w[IVX] *= tmols;
  w[IVY] *= tmols;
  w[IVZ] *= tmols;

  Real cvd = Constants::Rgas / (pthermo->GetGammadRef() - 1.);
  w[IEN] = tmols * (cvd * tem * fsig + LE);

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) w[n] *= tmols;

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) c[n] *= tmols;

  // tracer
  Real pd = xd / xgas * pres;
  Real rhod = pd / (tem * pthermo->GetRd());

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= rhod;

  SetType(Type::MoleConc);
}

void AirParcel::moleConcentrationToMoleFraction() {
  auto pthermo = Thermodynamics::GetInstance();
  Real tmols = w[IDN], fsig = w[IDN], xgas = 1., LE = 0.;

#pragma omp simd reduction(+ : tmols, fsig, LE)
  for (int n = 0; n < NCLOUD; ++n) {
    tmols += c[n];
    fsig += pthermo->GetCvRatioMole(n + 1 + NVAPOR) * c[n];
    LE += -pthermo->GetLatentEnergyMole(n + 1 + NVAPOR) * c[n];
  }

#pragma omp simd reduction(+ : tmols, fsig)
  for (int n = 1; n <= NVAPOR; ++n) {
    tmols += w[n];
    fsig += pthermo->GetCvRatioMole(n) * w[n];
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    c[n] /= tmols;
    xgas += -c[n];
  }

  Real xd = xgas;
#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    w[n] /= tmols;
    xd += -w[n];
  }

  w[IVX] /= tmols;
  w[IVY] /= tmols;
  w[IVZ] /= tmols;

  Real cvd = Constants::Rgas / (pthermo->GetGammadRef() - 1.);
  w[IDN] = (w[IEN] - LE) / (cvd * fsig);
  w[IPR] = xgas * tmols * Constants::Rgas * w[IDN];

  // tracer
  Real pd = xd / xgas * w[IPR];
  Real inv_rhod = (w[IDN] * pthermo->GetRd()) / pd;

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) x[n] *= inv_rhod;

  SetType(Type::MoleFrac);
}

void AirParcel::massFractionToMoleConcentration() {
  throw NotImplementedError("Variable::massFractionToMoleConcentration");
}

void AirParcel::moleConcentrationToMassFraction() {
  throw NotImplementedError("Variable::moleConcentrationToMassFraction");
}

void AirParcel::massConcentrationToMoleConcentration() {
  throw NotImplementedError("Variable::massConcentrationToMoleConcentration");
}

void AirParcel::moleConcentrationToMassConcentration() {
  throw NotImplementedError("Variable::moleConcentrationToMassConcentration");
}
