// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

// application
#include <application/exceptions.hpp>

// climath
#include <climath/core.h>

// harp
#include "harp/radiation.hpp"
#include "harp/radiation_band.hpp"

// inversion
#include "inversion/inversion.hpp"
#include "inversion/inversion_helper.hpp"

// snap
#include "snap/decomposition/decomposition.hpp"
#include "snap/implicit/implicit_solver.hpp"
#include "snap/thermodynamics/thermodynamics.hpp"

// canoe
#include "impl.hpp"
#include "index_map.hpp"

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // index map
  IndexMap::InitFromAthenaInput(pin);

  // thermodynamics
  Thermodynamics::InitFromAthenaInput(pin);

  // decomposition
  pdec = std::make_shared<Decomposition>(pmb);

  // implicit methods
  phevi = std::make_shared<ImplicitSolver>(pmb, pin);

  // cloud
  pcloud = std::make_shared<Cloud>(pmb, pin);

  // chemistry
  pchem = std::make_shared<Chemistry>(pmb, pin);

  // tracer
  ptracer = std::make_shared<Tracer>(pmb, pin);

  // radiation
  prad = std::make_shared<Radiation>(pmb, pin);

  // inversion queue
  fitq = Inversion::NewInversionQueue(pmb, pin);

  // reference pressure
#ifdef HYDROSTATIC
  reference_pressure_ = pin->GetReal("mesh", "ReferencePressure");

  // pressure scale height
  pressure_scale_height_ = pin->GetReal("mesh", "PressureScaleHeight");
#else
  reference_pressure_ = 1.0;
  pressure_scale_height_ = 1.0;
#endif  // HYDROSTATIC
}

MeshBlock::Impl::~Impl() {
  Thermodynamics::Destroy();
  IndexMap::Destroy();
}

void MeshBlock::Impl::GatherFromPrimitive(Variable *var, int k, int j,
                                          int i) const {
  auto &pmb = pmy_block_;
  auto mytype = var->GetType();
  auto phydro = pmb->phydro;
  var->SetType(Variable::Type::MassFrac);

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) var->w[n] = phydro->w(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) var->c[n] = pcloud->w(n, k, j, i);

  // scale mass fractions
  Real rho = var->w[IDN], xd = 1.;

#pragma omp simd reduction(+ : xd)
  for (int n = 1; n <= NVAPOR; ++n) xd += -var->w[n];
  Real rhod = var->w[IDN] * xd;

#pragma omp simd reduction(+ : rho)
  for (int n = 0; n < NCLOUD; ++n) rho += rhod * var->c[n];

  Real inv_rho = 1.0 / rho;

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n) {
    var->w[n] *= rhod * inv_rho;
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    var->c[n] *= rhod * inv_rho;
  }

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n)
    var->q[n] = pchem->w(n, k, j, i) * rhod * inv_rho;

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n)
    var->x[n] = ptracer->w(n, k, j, i) * rhod * inv_rho;

  var->w[IDN] = rho;

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::GatherFromConserved(Variable *var, int k, int j,
                                          int i) const {
  auto &pmb = pmy_block_;
  auto mytype = var->GetType();
  auto phydro = pmb->phydro;
  var->SetType(Variable::Type::MassConc);

  auto pthermo = Thermodynamics::GetInstance();

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) var->w[n] = phydro->u(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) var->c[n] = pcloud->u(n, k, j, i);

  Real rho = 0., cvt = 0.;
#pragma omp simd reduction(+ : rho, cvt)
  for (int n = 0; n <= NVAPOR; ++n) {
    rho += var->w[n];
    cvt += var->w[n] * pthermo->GetCvMassRef(n);
  }

  // TODO(cli): not correct for cubed sphere
  Real KE =
      0.5 / rho * (sqr(var->w[IVX]) + sqr(var->w[IVY]) + sqr(var->w[IVZ]));

  Real tem = (var->w[IEN] - KE) / cvt;

#pragma omp simd reduction(+ : cvt)
  for (int n = 0; n < NCLOUD; ++n) {
    cvt += var->c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR);
  }

  var->w[IEN] = cvt * tem;

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n) var->q[n] = pchem->u(n, k, j, i);

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) var->x[n] = ptracer->u(n, k, j, i);

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::DistributeToPrimitive(Variable const &var_in, int k,
                                            int j, int i) {
  Variable *var;
  auto &pmb = pmy_block_;
  auto phydro = pmb->phydro;

  if (var_in.GetType() != Variable::Type::MassFrac) {
    var = new Variable(var_in);
    var->ConvertToMassFraction();
  } else {
    var = const_cast<Variable *>(&var_in);
  }

  // scale mass fractions back
  Real rho = var->w[IDN];
  Real rhod = rho, rhog = rho;

#pragma omp simd reduction(+ : rhod, rhog)
  for (int n = 0; n < NCLOUD; ++n) {
    rhod += -rho * var->c[n];
    rhog += -rho * var->c[n];
  }

#pragma omp simd reduction(+ : rhod)
  for (int n = 1; n <= NVAPOR; ++n) rhod += -rho * var->w[n];

  Real inv_rhod = 1.0 / rhod;
  Real inv_rhog = 1.0 / rhog;

  phydro->w(IDN, k, j, i) = rhog;

#pragma omp simd
  for (int n = 1; n <= NVAPOR; ++n)
    phydro->w(n, k, j, i) = var->w[n] * rho * inv_rhog;

#pragma omp simd
  for (int n = 1 + NVAPOR; n < NHYDRO; ++n) phydro->w(n, k, j, i) = var->w[n];

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n)
    pcloud->w(n, k, j, i) = var->c[n] * rho * inv_rhod;

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n)
    pchem->w(n, k, j, i) = var->q[n] * rho * inv_rhod;

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n)
    ptracer->w(n, k, j, i) = var->x[n] * rho * inv_rhod;

  if (var_in.GetType() != Variable::Type::MassFrac) {
    delete var;
  }
}

void MeshBlock::Impl::DistributeToConserved(Variable const &var_in, int k,
                                            int j, int i) {
  Variable *var;
  auto &pmb = pmy_block_;
  auto phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();

  if (var_in.GetType() != Variable::Type::MassConc) {
    var = new Variable(var_in);
    var->ConvertToMassConcentration();
  } else {
    var = const_cast<Variable *>(&var_in);
  }

#pragma omp simd
  for (int n = 0; n < NHYDRO; ++n) phydro->u(n, k, j, i) = var->w[n];

  Real rhog = 0., cvt = 0.;

#pragma omp simd
  for (int n = 0; n <= NVAPOR; ++n) {
    rhog += var->w[n];
    cvt += var->w[n] * pthermo->GetCvMassRef(n);

    // TODO(cli): not correct for cubed sphere
    phydro->u(IEN, k, j, i) +=
        0.5 * var->w[n] * sqr(var->w[IVX] + var->w[IVY] + var->w[IVZ]);
  }

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    pcloud->u(n, k, j, i) = var->c[n];
    cvt += var->c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR);
  }

  Real tem = var->w[IEN] / cvt;

#pragma omp simd
  for (int n = 0; n < NCLOUD; ++n) {
    phydro->u(IEN, k, j, i) -=
        var->c[n] * pthermo->GetCvMassRef(n + 1 + NVAPOR) * tem;
  }

#pragma omp simd
  for (int n = 0; n < NCHEMISTRY; ++n) pchem->u(n, k, j, i) = var->q[n];

#pragma omp simd
  for (int n = 0; n < NTRACER; ++n) ptracer->u(n, k, j, i) = var->x[n];

  if (var_in.GetType() != Variable::Type::MassConc) {
    delete var;
  }
}

int find_pressure_level_lesser(Real pres, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie) {
  for (int i = is; i <= ie; ++i)
    if (w(IPR, k, j, i) < pres) return i;

  return ie + 1;
}
