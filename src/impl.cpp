// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

// application
#include <application/exceptions.hpp>

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

MeshBlock::Impl::Impl(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  du.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);

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

MeshBlock::Impl::~Impl() { Thermodynamics::Destroy(); }

void MeshBlock::Impl::GatherFromPrimitive(Variable *var, int k, int j, int i) {
  auto mytype = var->GetType();
  var->SetType(Variable::Type::MassFrac);

  for (int n = 0; n < NHYDRO; ++n)
    var->w[n] = pmy_block_->phydro->w(n, k, j, i);

  for (int n = 0; n < NCLOUD; ++n)
    var->c[n] = pmy_block_->pimpl->pcloud->w(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n)
    var->x[n] = pmy_block_->pimpl->ptracer->w(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n)
    var->q[n] = pmy_block_->pimpl->pchem->w(n, k, j, i);

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::GatherFromConserved(Variable *var, int k, int j, int i) {
  auto mytype = var->GetType();
  var->SetType(Variable::Type::MassConc);

  for (int n = 0; n < NHYDRO; ++n)
    var->w[n] = pmy_block_->phydro->u(n, k, j, i);

  for (int n = 0; n < NCLOUD; ++n)
    var->c[n] = pmy_block_->pimpl->pcloud->u(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n)
    var->x[n] = pmy_block_->pimpl->ptracer->u(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n)
    var->q[n] = pmy_block_->pimpl->pchem->u(n, k, j, i);

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::DistributeToPrimitive(Variable const& var_in, int k, int j, int i) {
  Variable *var;

  if (var_in.GetType() != Variable::Type::MassFrac) {
    var = new Variable(var_in);
    var->ConvertToMassFraction();
  } else {
    var = const_cast<Variable*>(&var_in);
  }

  for (int n = 0; n < NHYDRO; ++n) pmy_block_->phydro->w(n, k, j, i) = var->w[n];

  for (int n = 0; n < NCLOUD; ++n)
    pmy_block_->pimpl->pcloud->w(n, k, j, i) = var->c[n];

  for (int n = 0; n < NCHEMISTRY; ++n)
    pmy_block_->pimpl->pchem->w(n, k, j, i) = var->q[n];

  for (int n = 0; n < NTRACER; ++n)
    pmy_block_->pimpl->ptracer->w(n, k, j, i) = var->x[n];

  if (var_in.GetType() != Variable::Type::MassFrac) {
    delete var;
  }
}

void MeshBlock::Impl::DistributeToConserved(Variable const& var_in, int k, int j, int i) {
  Variable *var;

  if (var_in.GetType() != Variable::Type::MassConc) {
    var = new Variable(var_in);
    var->ConvertToMassConcentration();
  } else {
    var = const_cast<Variable*>(&var_in);
  }

  for (int n = 0; n < NHYDRO; ++n) pmy_block_->phydro->u(n, k, j, i) = var->w[n];

  for (int n = 0; n < NCLOUD; ++n)
    pmy_block_->pimpl->pcloud->u(n, k, j, i) = var->c[n];

  for (int n = 0; n < NCHEMISTRY; ++n)
    pmy_block_->pimpl->pchem->u(n, k, j, i) = var->q[n];

  for (int n = 0; n < NTRACER; ++n)
    pmy_block_->pimpl->ptracer->u(n, k, j, i) = var->x[n];

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
