// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.h>

// chemistry
#include "chemistry.hpp"

Chemistry::Chemistry(MeshBlock *pmb, ParameterInput *pin) {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("flask");
  app->Log("Initialize Chemistry");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD, NCHEMISTRY);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD, NCHEMISTRY);
}

Chemistry::~Chemistry() {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("flask");
  app->Log("Destroy Chemistry");
}
