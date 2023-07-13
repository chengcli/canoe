/** @file communicator.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Apr 14, 2022 11:48:11 EDT
 * @bug No known bugs.
 */

// C/C++
#include <iostream>

// athena
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// utils
#include "communicator.hpp"

Communicator::Communicator(MeshBlock *pmb) : pmy_block_(pmb) {
  pmb->pdebug->Enter("Communicator");
  color_ = new int[Globals::nranks];
  brank_ = new int[Globals::nranks];
#ifdef MPI_PARALLEL
  comm_ = MPI_COMM_WORLD;
#endif
  pmb->pdebug->Leave();
}

Communicator::~Communicator() {
  delete[] color_;
  delete[] brank_;

#ifdef MPI_PARALLEL
  if (comm_ != MPI_COMM_WORLD) MPI_Comm_free(&comm_);
#endif
}

int Communicator::getRank(CoordinateDirection dir) const {
  int r = 0;
  int b = brank_[Globals::my_rank];
  while (b != -1) {
    r++;
    b = brank_[b];
  }
  return r;
}

void Communicator::gatherData(Real *send, Real *recv, int size) const {
#ifdef MPI_PARALLEL
  MPI_Allgather(send, size, MPI_ATHENA_REAL, recv, size, MPI_ATHENA_REAL,
                comm_);
#else
  memcpy(recv, send, size * sizeof(Real));
#endif
}

//! \warning this function is unstable to some system.
void Communicator::gatherDataInPlace(Real *recv, int size) const {
#ifdef MPI_PARALLEL
  MPI_Allgather(MPI_IN_PLACE, 0, 0, recv, size, MPI_ATHENA_REAL, comm_);
#endif
}

NeighborBlock const *Communicator::findBotNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pbot = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == -1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0)) pbot = nb;
  }
  return pbot;
}

NeighborBlock const *Communicator::findTopNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *ptop = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0)) ptop = nb;
  }
  return ptop;
}

NeighborBlock const *Communicator::findLeftNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pleft = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 0) && (nb->ni.ox2 == -1) && (nb->ni.ox3 == 0))
      pleft = nb;
  }
  return pleft;
}

NeighborBlock const *Communicator::findRightNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pright = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 0))
      pright = nb;
  }
  return pright;
}

NeighborBlock const *Communicator::findBackNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pback = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 0) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == -1))
      pback = nb;
  }
  return pback;
}

NeighborBlock const *Communicator::findFrontNeighbor() const {
  MeshBlock *pmb = pmy_block_;
  NeighborBlock *pfront = nullptr;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 1))
      pfront = nb;
  }
  return pfront;
}

void Communicator::setColor(CoordinateDirection dir) {
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
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_, 1, MPI_INT,
                MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  std::fill(color_, color_ + Globals::nranks, -1);
  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    // color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color_[i] == -1) color_[i] = c++;
    } else
      color_[i] = color_[brank_[i]];
  }

#ifdef MPI_PARALLEL
  if (comm_ != MPI_COMM_WORLD) MPI_Comm_free(&comm_);
  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &comm_);
#endif
}
