// C/C++ headers
#include <cstring>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "diagnostics.hpp"

Diagnostics::Diagnostics(MeshBlock *pmb, ParameterInput *pin)
    : myname("HEAD"),
      type(""),
      grid(""),
      varname("HEAD"),
      long_name(""),
      units(""),
      prev(nullptr),
      next(nullptr),
      ncycle(0),
      pmy_block_(pmb) {
  pmb->pdebug->Enter("Diagnostics List");
  std::stringstream msg;
  char cstr[80];
  std::string diag_names = pin->GetOrAddString("problem", "diagnostics", "");
  std::strcpy(cstr, diag_names.c_str());
  char *p = std::strtok(cstr, " ,");
  while (p != NULL) {
    std::string name(p);
    if (name == "div") {  // 1.
      AddDiagnostics(Divergence(pmb));
    } else if (name == "curl")  // 2.
      AddDiagnostics(Curl(pmb));
    else if (name == "mean")  // 3.
      AddDiagnostics(HydroMean(pmb));
    else if (name == "tempa")  // 4.
      AddDiagnostics(TemperatureAnomaly(pmb));
    else if (name == "presa")  // 5.
      AddDiagnostics(PressureAnomaly(pmb));
    else if (name == "eddyflux")  // 6.
      AddDiagnostics(EddyFlux(pmb));
    else if (name == "hydroflux")  // 7.
      AddDiagnostics(HydroFlux(pmb));
    else if (name == "div_h")  // 8.
      AddDiagnostics(HorizontalDivergence(pmb));
    else if (name == "b")  // 9.
      AddDiagnostics(Buoyancy(pmb));
    else if (name == "radflux")  // 10.
      AddDiagnostics(RadiativeFlux(pmb));
    else if (name == "am")  // 11.
      AddDiagnostics(AngularMomentum(pmb));
    else if (name == "eke")  // 12.
      AddDiagnostics(EddyKineticEnergy(pmb));
    else if (name == "tendency")  // 13.
      AddDiagnostics(Tendency(pmb));
    else {
      msg << "### FATAL ERROR in function Diagnostics::Diagnostics" << std::endl
          << "Diagnostic variable " << name << " not defined";
      ATHENA_ERROR(msg);
    }
    // msg << "- add diagnostics " + name << std::endl;
    p = std::strtok(NULL, " ,");
  }

#if MAIN_TASKLIST != InversionTaskList
  if (NGHOST < 2) {
    msg << "### FATAL ERROR in function Diagnostics::Diagnostics" << std::endl
        << "Most diagnostic variables require at least 2 ghost cells";
    ATHENA_ERROR(msg);
  }
#endif

  // pmb->pdebug->WriteMessage(msg.str());
  pmb->pdebug->Leave();
}

Diagnostics::Diagnostics(MeshBlock *pmb, std::string name)
    : myname(name),
      varname(name),
      prev(nullptr),
      next(nullptr),
      ncycle(0),
      pmy_block_(pmb) {
  pmb->pdebug->Enter("Diagnostics-" + name);
  std::stringstream msg;
  ncells1_ = pmb->block_size.nx1 + 2 * (NGHOST);
  ncells2_ = 1;
  ncells3_ = 1;
  if (pmb->pmy_mesh->f2) ncells2_ = pmb->block_size.nx2 + 2 * (NGHOST);
  if (pmb->pmy_mesh->f3) ncells3_ = pmb->block_size.nx3 + 2 * (NGHOST);

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
  total_vol_.NewAthenaArray(ncells1_);
  brank_.resize(Globals::nranks);
  color_.resize(Globals::nranks);

  for (int i = 0; i < Globals::nranks; ++i) {
    brank_[i] = -1;
    color_[i] = -1;
  }

  pmb->pdebug->Leave();
}

Diagnostics::~Diagnostics() {
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
}

Diagnostics *Diagnostics::operator[](std::string name) {
  Diagnostics *p = this;
  while (p != nullptr) {
    if (p->myname == name) return p;
    p = p->next;
  }
  return p;
}

void Diagnostics::setColor_(int *color, CoordinateDirection dir) {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock bblock, tblock;
  pmb->FindNeighbors(dir, bblock, tblock);

  if (dir == X1DIR) {
    if (pmb->block_size.x1min <= pmb->pmy_mesh->mesh_size.x1min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x1max >= pmb->pmy_mesh->mesh_size.x1max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else if (dir == X2DIR) {
    if (pmb->block_size.x2min <= pmb->pmy_mesh->mesh_size.x2min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x2max >= pmb->pmy_mesh->mesh_size.x2max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  } else {  // X3DIR
    if (pmb->block_size.x3min <= pmb->pmy_mesh->mesh_size.x3min) {
      bblock.snb.gid = -1;
      bblock.snb.rank = -1;
    }
    if (pmb->block_size.x3max >= pmb->pmy_mesh->mesh_size.x3max) {
      tblock.snb.gid = -1;
      tblock.snb.rank = -1;
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    // color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color[i] == -1) color[i] = c++;
    } else
      color[i] = color[brank_[i]];
  }
}

void Diagnostics::gatherAllData23_(AthenaArray<Real> &total_vol,
                                   AthenaArray<Real> &total_data) {
  MeshBlock *pmb = pmy_block_;

  // calculate total volume
  total_vol.ZeroClear();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
      for (int i = pmb->is; i <= pmb->ie; ++i) total_vol(i) += vol_(i);
    }

    // sum over all ranks
#ifdef MPI_PARALLEL
  MPI_Comm comm;
  std::fill(color_.data(), color_.data() + Globals::nranks, -1);
  setColor_(color_.data(), X2DIR);
  setColor_(color_.data(), X3DIR);
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &comm);
  // int size;
  // MPI_Comm_size(comm, &size);
  // std::cout << size << std::endl;
  MPI_Allreduce(MPI_IN_PLACE, total_vol.data(), total_vol.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, total_data.data(), total_data.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, comm);
  MPI_Comm_free(&comm);
#endif
}
