// athena
#include <athena/bvals/bvals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <common.hpp>

namespace NeighborExchangerHelper {

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

NeighborBlock const *find_right_neighbor(MeshBlock const *pmb) const {
  NeighborBlock *pright = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 0))
      pright = nb;
  }

  return pright;
}

NeighborBlock const *find_back_neighbor(MeshBlock const *pmb) const {
  NeighborBlock *pback = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 0) && (nb->ni.ox2 == 0) && (nb->ni.ox3 == -1))
      pback = nb;
  }

  return pback;
}

NeighborBlock const *find_front_neighbor(MeshBlock const *pmb) const {
  NeighborBlock *pfront = nullptr;

  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock *nb = pmb->pbval->neighbor + n;
    if ((nb->ni.ox1 == 1) && (nb->ni.ox2 == 1) && (nb->ni.ox3 == 1))
      pfront = nb;
  }

  return pfront;
}

void find_neighbors(MeshBlock const *pmb, CoordinateID dir,
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

}  // namespace NeighborExchangerHelper
