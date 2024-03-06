// external
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <virtual_groups.hpp>

// scm
#include "single_column.hpp"

SingleColumn::SingleColumn(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  // auto app = Application::GetInstance();
  // app->InstallMonitor("single_column", "single_column.out",
  // "single_column.out"); auto log = app->GetMonitor("single_column");

  Application::Logger app("single_column");
  app->Log("Initialize SingleColumn");

  SetPar("den_tol",
         pin->GetOrAddReal("convective_adjustment", "den_tol", 0.001));
  SetPar("rel_tol",
         pin->GetOrAddReal("convective_adjustment", "rel_tol", 1.0e-4));
  SetPar("max_iter",
         pin->GetOrAddInteger("convective_adjustment", "max_iter", 10));

  // allocate scratch arrays
  vol_.NewAthenaArray(pmb->ncells1);
}

SingleColumn::~SingleColumn() {
  Application::Logger app("single_column");
  app->Log("Destroy SingleColumn");
}
