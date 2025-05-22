// C/C++ headers
#include <algorithm>
#include <string>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.h>

// snap
#include "turbulence_model.hpp"

// constructor, initializes data structures and parameters

TurbulenceModel::TurbulenceModel(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block(pmb) {
  if (NTURBULENCE == 0) return;

  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD + NCHEMISTRY + NTRACER,
                         NTURBULENCE);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD + NCHEMISTRY + NTRACER,
                         NTURBULENCE);

  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  mut.NewAthenaArray(ncells3, ncells2, ncells1);
}

TurbulenceModel::~TurbulenceModel() {
  if (NTURBULENCE == 0) return;

  Application::Logger app("snap");
  app->Log("Destroy TurbulenceModel");
}

TurbulenceModelPtr TurbulenceFactory::Create(MeshBlock *pmb,
                                             ParameterInput *pin) {
  Application::Logger app("snap");
  app->Log("Create turbulenceModel");

  TurbulenceModelPtr pturb = nullptr;

  if (pin->DoesParameterExist("hydro", "turbulence")) {
    std::string turbulence_model = pin->GetString("hydro", "turbulence");
    if (turbulence_model == "none") {
      pturb = std::make_shared<TurbulenceModel>(pmb, pin);
    } else if (turbulence_model == "kepsilon") {
      pturb = std::make_shared<KEpsilonTurbulence>(pmb, pin);
      if (NTURBULENCE < 2)
        throw NotImplementedError(
            "NTURBULENCE must be at least 2 for k-epsilon model");
    } else {
      throw NotImplementedError(turbulence_model);
    }
  }

  return pturb;
}
