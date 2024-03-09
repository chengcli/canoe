#ifndef SRC_EXCHANGER_LINEAR_EXCHANGER_HPP_
#define SRC_EXCHANGER_LINEAR_EXCHANGER_HPP_

// athena
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// exchanger
#include "exchanger.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

template <typename T>
LinearExchanger<T>::LinearExchanger() {
#ifdef MPI_PARALLEL
  color_.resize(Globals::nranks);
  brank_.resize(Globals::nranks);
#else
  color_.resize(1);
  brank_.resize(1);
  color_[0] = 0;
  brank_[0] = -1;
#endif  // MPI_PARALLEL
}

template <typename T>
int LinearExchanger<T>::GetRankInGroup() const {
  int r = 0;
  int b = brank_[Globals::my_rank];
  while (b != -1) {
    r++;
    b = brank_[b];
  }
  return r;
}

template <typename T>
void LinearExchanger<T>::Regroup(MeshBlock const *pmb,
                                 CoordinateDirection dir) {
#ifdef MPI_PARALLEL
  NeighborBlock bblock, tblock;
  ExchangerHelper::find_neighbors(pmb, dir, &bblock, &tblock);

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

  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  std::fill(color_.begin(), color_.end(), -1);
  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    // color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color_[i] == -1) color_[i] = c++;
    } else
      color_[i] = color_[brank_[i]];
  }

#ifdef MPI_PARALLEL
  if (mpi_comm_ != MPI_COMM_WORLD) {
    MPI_Comm_free(&mpi_comm_);
  }

  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &mpi_comm_);
#endif
}

#endif  // SRC_EXCHANGER_LINEAR_EXCHANGER_HPP_
