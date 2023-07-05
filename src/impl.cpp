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

  // implicit methodsphydro
  phevi = std::make_shared<ImplicitSolver>(pmb, pin);

  // tracer
  ptracer = std::make_shared<Tracer>(pmb, pin);

  // chemistry
  pchem = std::make_shared<Chemistry>(pmb, pin);

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

void MeshBlock::Impl::GatherPrimitive(Variable *var, int k, int j, int i) {
  for (int n = 0; n < NHYDRO; ++n)
    var->w[n] = pmy_block_->phydro->w(n, k, j, i);

  // for (int n = 0; n < NCLOUD; ++n)
  //   var->c[n] = pmy_block_->pimpl->pcloud->w(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n)
    var->x[n] = pmy_block_->pimpl->ptracer->w(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n)
    var->q[n] = pmy_block_->pimpl->pchem->w(n, k, j, i);

  var->SetType(Variable::Type::MassFrac);
}

void MeshBlock::Impl::GatherConserved(Variable *var, int k, int j, int i) {
  for (int n = 0; n < NHYDRO; ++n)
    var->w[n] = pmy_block_->phydro->u(n, k, j, i);

  // for (int n = 0; n < NCLOUD; ++n)
  //   var->c[n] = pmy_block_->pimpl->pcloud->u(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n)
    var->x[n] = pmy_block_->pimpl->ptracer->u(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n)
    var->q[n] = pmy_block_->pimpl->pchem->u(n, k, j, i);

  var->SetType(Variable::Type::MassConc);
}

void MeshBlock::Impl::GatherMoleFraction(Variable *var, int k, int j, int i) {}

void MeshBlock::Impl::DistributePrimitive(Variable &var, int k, int j, int i) {
  var.ConvertToMassFraction();

  for (int n = 0; n < NHYDRO; ++n) pmy_block_->phydro->w(n, k, j, i) = var.w[n];

  // for (int n = 0; n < NCLOUD; ++n)
  //   pmy_block_->pimpl->pcloud->w(n, k, j, i) = var.c[n];

  for (int n = 0; n < NTRACER; ++n)
    pmy_block_->pimpl->ptracer->w(n, k, j, i) = var.x[n];

  for (int n = 0; n < NCHEMISTRY; ++n)
    pmy_block_->pimpl->pchem->w(n, k, j, i) = var.q[n];
}

void MeshBlock::Impl::DistributeConserved(Variable &var, int k, int j, int i) {
  var.ConvertToMassConcentration();

  for (int n = 0; n < NHYDRO; ++n) pmy_block_->phydro->u(n, k, j, i) = var.w[n];

  // for (int n = 0; n < NCLOUD; ++n)
  //   pmy_block_->pimpl->pcloud->u(n, k, j, i) = var.c[n];

  for (int n = 0; n < NTRACER; ++n)
    pmy_block_->pimpl->ptracer->u(n, k, j, i) = var.x[n];

  for (int n = 0; n < NCHEMISTRY; ++n)
    pmy_block_->pimpl->pchem->u(n, k, j, i) = var.q[n];
}

int find_pressure_level_lesser(Real pres, AthenaArray<Real> const &w, int k,
                               int j, int is, int ie) {
  for (int i = is; i <= ie; ++i)
    if (w(IPR, k, j, i) < pres) return i;

  return ie + 1;
}
