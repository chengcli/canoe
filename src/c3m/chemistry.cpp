// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// chemistry
#include "chemistry.hpp"

Kinetics1D::Kinetics1D(MeshBlock *pmb, ParameterInput *pin) : pmy_block_(pmb) {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("c3m");
  app->Log("Initialize Kinetics1D");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD, NCHEMISTRY);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD, NCHEMISTRY);
}

void Kinetics1D::fromMeshBlock(MeshBlock *pmb, int k, int j) {
  if (NCHEMISTRY == 0) return;

  resize(nComponents(), pmb->block_size.nx1);

  pmy_block_ = pmb;
  w.InitWithShallowSlice(pmb->pscalars->r, 4, NCLOUD, NCHEMISTRY);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, NCLOUD, NCHEMISTRY);
}

Kinetics1D::~Kinetics1D() {
  if (NCHEMISTRY == 0) return;

  Application::Logger app("c3m");
  app->Log("Destroy Kinetics1D");
}
