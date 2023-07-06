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
  var->SetType(Variable::Type::MassFrac);

  for (int n = 0; n < NHYDRO; ++n) var->w[n] = pmb->phydro->w(n, k, j, i);

  for (int n = 0; n < NCLOUD; ++n) var->c[n] = pcloud->w(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n) var->x[n] = ptracer->w(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n) var->q[n] = pchem->w(n, k, j, i);

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::GatherFromConserved(Variable *var, int k, int j,
                                          int i) const {
  auto &pmb = pmy_block_;
  auto mytype = var->GetType();
  var->SetType(Variable::Type::MassConc);

  for (int n = 0; n < NHYDRO; ++n) var->w[n] = pmb->phydro->u(n, k, j, i);

  for (int n = 0; n < NCLOUD; ++n) var->c[n] = pcloud->u(n, k, j, i);

  for (int n = 0; n < NTRACER; ++n) var->x[n] = ptracer->u(n, k, j, i);

  for (int n = 0; n < NCHEMISTRY; ++n) var->q[n] = pchem->u(n, k, j, i);

  var->ConvertTo(mytype);
}

void MeshBlock::Impl::DistributeToPrimitive(Variable const &var_in, int k,
                                            int j, int i) {
  Variable *var;
  auto &pmb = pmy_block_;

  if (var_in.GetType() != Variable::Type::MassFrac) {
    var = new Variable(var_in);
    var->ConvertToMassFraction();
  } else {
    var = const_cast<Variable *>(&var_in);
  }

  for (int n = 0; n < NHYDRO; ++n) pmb->phydro->w(n, k, j, i) = var->w[n];

  for (int n = 0; n < NCLOUD; ++n) pcloud->w(n, k, j, i) = var->c[n];

  for (int n = 0; n < NCHEMISTRY; ++n) pchem->w(n, k, j, i) = var->q[n];

  for (int n = 0; n < NTRACER; ++n) ptracer->w(n, k, j, i) = var->x[n];

  if (var_in.GetType() != Variable::Type::MassFrac) {
    delete var;
  }
}

void MeshBlock::Impl::DistributeToConserved(Variable const &var_in, int k,
                                            int j, int i) {
  Variable *var;
  auto &pmb = pmy_block_;

  if (var_in.GetType() != Variable::Type::MassConc) {
    var = new Variable(var_in);
    var->ConvertToMassConcentration();
  } else {
    var = const_cast<Variable *>(&var_in);
  }

  for (int n = 0; n < NHYDRO; ++n) pmb->phydro->u(n, k, j, i) = var->w[n];

  for (int n = 0; n < NCLOUD; ++n) pcloud->u(n, k, j, i) = var->c[n];

  for (int n = 0; n < NCHEMISTRY; ++n) pchem->u(n, k, j, i) = var->q[n];

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
