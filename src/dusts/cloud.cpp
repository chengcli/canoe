// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// dusts
#include "cloud.hpp"

Cloud::Cloud(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NCLOUD == 0) return;

  Application::Logger app("dusts");
  app->Log("Initialize Cloud");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, 0, NCLOUD);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, 0, NCLOUD);
}

Cloud::~Cloud() {
  if (NCLOUD == 0) return;

  Application::Logger app("dusts");
  app->Log("Destroy Cloud");
}
