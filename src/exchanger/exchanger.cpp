// athena
#include <athena/bvals/bvals.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include "exchanger.hpp"

ExchangerBase::ExchangerBase() {
#ifdef MPI_PARALLEL
  mpi_comm_ = MPI_COMM_WORLD;
  color_.resize(Globals::nranks);
  brank_.resize(Globals::nranks);

  for (int i = 0; i < Globals::nranks; ++i) {
    color_[i] = -1;
    brank_[i] = -1;
  }
#else
  color_.resize(1);
  brank_.resize(1);
  color_[0] = 0;
  brank_[0] = -1;
#endif  // MPI_PARALLEL
}

void ExchangerBase::setColor(MeshBlock const *pmb, CoordinateDirection dir) {
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

  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    // color[i] = brank_[i] == -1 ? color[i] : color[brank_[i]];
    if (brank_[i] == -1) {
      if (color_[i] == -1) color_[i] = c++;
    } else {
      color_[i] = color_[brank_[i]];
    }
  }

#ifdef MPI_PARALLEL
  if (mpi_comm_ != MPI_COMM_WORLD) {
    MPI_Comm_free(&mpi_comm_);
  }

  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank], Globals::my_rank,
                 &mpi_comm_);
#endif
}

int ExchangerBase::GetRankInGroup() const {
  int r = 0;
  int b = brank_[Globals::my_rank];
  while (b != -1) {
    r++;
    b = brank_[b];
  }
  return r;
}

namespace ExchangerHelper {

int mpi_tag_ub;

int create_mpi_tag(int lid, int tid, std::string name) {
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, 0);

  std::string str = name + std::to_string(tag);
  return std::hash<std::string>{}(str) % (mpi_tag_ub);
}

NeighborBlock const *find_bot_neighbor(MeshBlock const *pmb) {
  NeighborBlock *pbot = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == -1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0)) pbot = nb;
  }

  return pbot;
}

NeighborBlock const *find_top_neighbor(MeshBlock const *pmb) {
  NeighborBlock *ptop = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == 0)) ptop = nb;
  }

  return ptop;
}

NeighborBlock const *find_left_neighbor(MeshBlock const *pmb) {
  NeighborBlock *pleft = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 0) && (nb->ni.ox2 == -1) && (nb->ni.ox3 == 0))
      pleft = nb;
  }

  return pleft;
}

NeighborBlock const *find_right_neighbor(MeshBlock const *pmb) {
  NeighborBlock *pright = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 0))
      pright = nb;
  }

  return pright;
}

NeighborBlock const *find_back_neighbor(MeshBlock const *pmb) {
  NeighborBlock *pback = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 0) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == -1))
      pback = nb;
  }

  return pback;
}

NeighborBlock const *find_front_neighbor(MeshBlock const *pmb) {
  NeighborBlock *pfront = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 1))
      pfront = nb;
  }

  return pfront;
}

void find_neighbors(MeshBlock const *pmb, CoordinateDirection dir,
                    NeighborBlock *bblock, NeighborBlock *tblock) {
  // set void bblock
  bblock->snb.gid = -1;
  bblock->snb.rank = -1;

  // set void tblock
  tblock->snb.gid = -1;
  tblock->snb.rank = -1;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (dir == X1DIR) {
      if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
        *bblock = nb;
      if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
        *tblock = nb;
    } else if (dir == X2DIR) {
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == -1) && (nb.ni.ox3 == 0))
        *bblock = nb;
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 1) && (nb.ni.ox3 == 0))
        *tblock = nb;
    } else {  // X3DIR
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == -1))
        *bblock = nb;
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 1))
        *tblock = nb;
    }
  }
}

}  // namespace ExchangerHelper
