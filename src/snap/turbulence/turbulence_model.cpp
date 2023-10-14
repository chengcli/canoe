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

// canoe
#include <configure.hpp>

// snap
#include "turbulence_model.hpp"

// constructor, initializes data structures and parameters

TurbulenceModel::TurbulenceModel(MeshBlock *pmb, ParameterInput *pin, int nvar)
    : pmy_block(pmb) {
  if (NTURBULENCE == 0) return;

  Application::Logger app("snap");
  app->Log("Initialize Turbulence");

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
  app->Log("Destroy Turbulence");
}
