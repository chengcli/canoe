// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// flask
#include "flask.hpp"

Flask::Flask(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("flask");
  app->Log("Initialize Flask");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD, NCHEMISTRY);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD, NCHEMISTRY);
}

Flask::~Flask() {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("flask");
  app->Log("Destroy Flask");
}
