// C/C++ headers
#include <cstring>
#include <sstream>
#include <stdexcept>

// application
#include <application/application.hpp>

// Athena++ headers
#include <athena/mesh/mesh.hpp>

// diagnostics
#include "diagnostics.hpp"

Diagnostics::Diagnostics(MeshBlock *pmb, std::string name) : NamedGroup(name) {
  Application::Logger app("main");
  app->Log("Initialize Diagnostics");

  ncells1_ = pmb->block_size.nx1 + 2 * (NGHOST);
  ncells2_ = 1;
  ncells3_ = 1;

  if (pmb->pmy_mesh->f2) {
    ncells2_ = pmb->block_size.nx2 + 2 * (NGHOST);
  }

  if (pmb->pmy_mesh->f3) {
    ncells3_ = pmb->block_size.nx3 + 2 * (NGHOST);
  }

  x1edge_.NewAthenaArray(ncells1_ + 1);
  x1edge_p1_.NewAthenaArray(ncells1_);
  x2edge_.NewAthenaArray(ncells1_ + 1);
  x2edge_p1_.NewAthenaArray(ncells1_);
  x3edge_.NewAthenaArray(ncells1_ + 1);
  x3edge_p1_.NewAthenaArray(ncells1_);

  x1area_.NewAthenaArray(ncells1_ + 1);
  x2area_.NewAthenaArray(ncells1_);
  x2area_p1_.NewAthenaArray(ncells1_);
  x3area_.NewAthenaArray(ncells1_);
  x3area_p1_.NewAthenaArray(ncells1_);

  vol_.NewAthenaArray(ncells1_);
  total_vol_.resize(ncells1_);
  total_area_.resize(ncells1_);
}

Diagnostics::~Diagnostics() {
  Application::Logger app("main");
  app->Log("Destroy Diagnostics");
}
