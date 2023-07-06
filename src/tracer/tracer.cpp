// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// tracer
#include "tracer.hpp"

Tracer::Tracer(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NTRACER == 0) return;

  Application::Logger app("tracer");
  app->Log("Initialize Tracer");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD + NCHEMISTRY, NTRACER);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD + NCHEMISTRY, NTRACER);
}

Tracer::~Tracer() {
  if (NTRACER == 0) return;

  Application::Logger app("tracer");
  app->Log("Destroy Tracer");
}
