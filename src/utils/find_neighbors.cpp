// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

void find_neighbors(MeshBlock *pmb, CoordinateDirection dir,
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
