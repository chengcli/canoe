// application
#include <application/exceptions.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/core.h>  // sqr

// canoe
#include "air_parcel.hpp"
#include "impl.hpp"

// microphysics
#include "microphysics/microphysics.hpp"

// chemistry
#include "c3m/chemistry.hpp"

// tracer
#include "tracer/tracer.hpp"

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

  Real fsig = w[IDN], feps = w[IDN];
  // vapors
#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n] * pthermo->GetCvRatioMass(n);
    feps += w[n] * pthermo->GetInvMuRatio(n);
  }

  // clouds
#pragma omp simd reduction(+ : fsig)
  for (int n = 0; n < NCLOUD; ++n) {
    fsig += c[n] * pthermo->GetCvRatioMass(n + 1 + NVAPOR);
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

#pragma omp simd reduction(+ : fsig, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    fsig += c[n] * (pthermo->GetCvRatioMass(n + 1 + NVAPOR) - 1.);
    feps += -c[n];
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

// thermodynamic functions
Real AirParcel::gammad() const { return 1.4; }

Real AirParcel::chi() const {
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel* air;

  if (mytype_ != AirParcel::Type::MoleFrac) {
    air = new AirParcel(*this);
    air->ToMoleFraction();
  } else {
    air = const_cast<AirParcel*>(this);
  }

  Real gammad = air->gammad();

  Real qsig = 1., feps = 1.;
#pragma omp simd reduction(+ : qsig)
  for (int n = 1; n <= NVAPOR; ++n) {
    qsig += air->w[n] * (pthermo->GetCpRatioMole(n) - 1.);
  }

#pragma omp simd reduction(+ : qsig, feps)
  for (int n = 0; n < NCLOUD; ++n) {
    feps += -air->c[n];
    qsig += air->c[n] * (pthermo->GetCpRatioMole(n + 1 + NVAPOR) - 1.);
  }

  if (mytype_ != AirParcel::Type::MoleFrac) {
    delete air;
  }

  return (gammad - 1.) / gammad / qsig;
}

Real AirParcel::theta(Real p0) const {
  AirParcel* air;

  if (mytype_ != AirParcel::Type::MoleFrac) {
    air = new AirParcel(*this);
    air->ToMoleFraction();
  } else {
    air = const_cast<AirParcel*>(this);
  }

  if (mytype_ != AirParcel::Type::MoleFrac) {
    delete air;
  }
  return air->w[IDN] * pow(p0 / air->w[IPR], air->chi());
}

namespace AirParcelHelper {

AirParcel gather_from_primitive(MeshBlock const* pmb, int k, int j, int i) {
  AirParcel air(AirParcel::Type::MassFrac);

  auto phydro = pmb->phydro;
  auto ptracer = pmb->pimpl->ptracer;
  auto pmicro = pmb->pimpl->pmicro;
  auto pchem = pmb->pimpl->pchem;

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) air.w[n] = phydro->w(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) air.c[n] = pmicro->w(n, k, j, i);

  // scale mass fractions
  Real rho = air.w[IDN], xd = 1.;

#pragma omp simd reduction(+ : xd)
  for (int n = 1; n <= NVAPOR; ++n) xd += -air.w[n];
  Real rhod = air.w[IDN] * xd;

#pragma omp simd reduction(+ : rho)
  for (int n = 0; n < NCLOUD; ++n) rho += rhod * air.c[n];

  Real inv_rho = 1.0 / rho;

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    air.w[n] *= air.w[IDN] * inv_rho;
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    air.c[n] *= rhod * inv_rho;
  }

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n)
    air.q[n] = pchem->w(n, k, j, i) * rhod * inv_rho;

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) air.x[n] = ptracer->w(n, k, j, i);

  air.w[IDN] = rho;

  return air;
}

AirParcel gather_from_conserved(MeshBlock const* pmb, int k, int j, int i) {
  AirParcel air(AirParcel::Type::MassConc);

  auto phydro = pmb->phydro;
  auto ptracer = pmb->pimpl->ptracer;
  auto pmicro = pmb->pimpl->pmicro;
  auto pchem = pmb->pimpl->pchem;

  auto pthermo = Thermodynamics::GetInstance();

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) air.w[n] = phydro->u(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) air.c[n] = pmicro->u(n, k, j, i);

  Real rho = 0., cvt = 0.;
#pragma omp simd reduction(+ : rho, cvt)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += air.w[n];
    cvt += air.w[n] * pthermo->GetCvMassRef(n);
  }

  //! \todo not correct for cubed sphere
  Real KE = 0.5 / rho * (sqr(air.w[IVX]) + sqr(air.w[IVY]) + sqr(air.w[IVZ]));

  Real tem = (air.w[IEN] - KE) / cvt;
  Real vx = air.w[IVX] / rho;
  Real vy = air.w[IVY] / rho;
  Real vz = air.w[IVZ] / rho;

#pragma omp simd reduction(+ : cvt)
  for (int n = 0; n < NCLOUD; ++n) {
    cvt += air.c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR);
    air.w[IVX] += air.c[n] * vx;
    air.w[IVY] += air.c[n] * vy;
    air.w[IVZ] += air.c[n] * vz;
  }

  air.w[IEN] = cvt * tem;

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n) air.q[n] = pchem->u(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) air.x[n] = ptracer->u(n, k, j, i);

  return air;
}

void distribute_to_primitive(MeshBlock* pmb, int k, int j, int i,
                             AirParcel const& air_in) {
  AirParcel* air;

  auto phydro = pmb->phydro;
  auto ptracer = pmb->pimpl->ptracer;
  auto pmicro = pmb->pimpl->pmicro;
  auto pchem = pmb->pimpl->pchem;

  if (air_in.GetType() != AirParcel::Type::MassFrac) {
    air = new AirParcel(air_in);
    air->ToMassFraction();
  } else {
    air = const_cast<AirParcel*>(&air_in);
  }

  // scale mass fractions back
  Real rho = air->w[IDN];
  Real rhod = rho, rhog = rho;

#pragma omp simd reduction(+ : rhod, rhog)
  for (int n = 0; n < NCLOUD; ++n) {
    rhod += -rho * air->c[n];
    rhog += -rho * air->c[n];
  }

#pragma omp simd reduction(+ : rhod)
  for (int n = 1; n <= NVAPOR; ++n) rhod += -rho * air->w[n];

  Real inv_rhod = 1.0 / rhod;
  Real inv_rhog = 1.0 / rhog;

  phydro->w(IDN, k, j, i) = rhog;

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n)
    phydro->w(n, k, j, i) = air->w[n] * rho * inv_rhog;

#pragma omp simd
  for (int n = 1 + NVAPOR; n < NHYDRO; ++n) phydro->w(n, k, j, i) = air->w[n];

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n)
    pmicro->w(n, k, j, i) = air->c[n] * rho * inv_rhod;

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n)
    pchem->w(n, k, j, i) = air->q[n] * rho * inv_rhod;

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) ptracer->w(n, k, j, i) = air->x[n];

  if (air_in.GetType() != AirParcel::Type::MassFrac) {
    delete air;
  }
}

void distribute_to_conserved(MeshBlock* pmb, int k, int j, int i,
                             AirParcel const& air_in) {
  AirParcel* air;

  auto phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();
  auto ptracer = pmb->pimpl->ptracer;
  auto pmicro = pmb->pimpl->pmicro;
  auto pchem = pmb->pimpl->pchem;

  if (air_in.GetType() != AirParcel::Type::MassConc) {
    air = new AirParcel(air_in);
    air->ToMassConcentration();
  } else {
    air = const_cast<AirParcel*>(&air_in);
  }

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) phydro->u(n, k, j, i) = air->w[n];

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) pmicro->u(n, k, j, i) = air->c[n];

  Real rho = 0., cvt = 0.;

#pragma omp simd reduction(+ : rho, cvt)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += air->w[n];
    cvt += air->w[n] * pthermo->GetCvMassRef(n);
  }

#pragma omp simd reduction(+ : rho, cvt)
  for (int n = 0; n < NCLOUD; ++n) {
    rho += air->c[n];
    cvt += air->c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR);
  }

  Real tem = air->w[IEN] / cvt;
  Real vx = air->w[IVX] / rho;
  Real vy = air->w[IVY] / rho;
  Real vz = air->w[IVZ] / rho;

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    phydro->u(IEN, k, j, i) -=
        air->c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR) * tem;
    phydro->u(IVX, k, j, i) -= air->c[n] * vx;
    phydro->u(IVY, k, j, i) -= air->c[n] * vy;
    phydro->u(IVZ, k, j, i) -= air->c[n] * vz;
  }

  for (int n = 0; n <= NVAPOR; ++n) {
    // TODO(cli): not correct for cubed sphere
    phydro->u(IEN, k, j, i) += 0.5 * air->w[n] * (sqr(vx) + sqr(vy) + sqr(vz));
  }

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n) pchem->u(n, k, j, i) = air->q[n];

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) ptracer->u(n, k, j, i) = air->x[n];

  if (air_in.GetType() != AirParcel::Type::MassConc) {
    delete air;
  }
}

AirColumn gather_from_primitive(MeshBlock const* pmb, int k, int j) {
  return gather_from_primitive(pmb, k, j, 0, pmb->ncells1 - 1);
}

AirColumn gather_from_conserved(MeshBlock const* pmb, int k, int j) {
  return gather_from_conserved(pmb, k, j, 0, pmb->ncells1 - 1);
}

void distribute_to_primitive(MeshBlock* pmb, int k, int j,
                             AirColumn const& ac) {
  distribute_to_primitive(pmb, k, j, 0, pmb->ncells1 - 1, ac);
}

void distribute_to_conserved(MeshBlock* pmb, int k, int j,
                             AirColumn const& ac) {
  distribute_to_conserved(pmb, k, j, 0, pmb->ncells1 - 1, ac);
}

}  // namespace AirParcelHelper
