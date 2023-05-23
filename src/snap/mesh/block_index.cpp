// Athena++ headers
#include <mesh/mesh.hpp>

// canoe headers
#include "block_index.hpp"

// one additional grid for face variables
BlockIndex::BlockIndex(int nk, int nj, int ni, int mygid, int mylid)
    : is(nghost),
      ie(nghost + ni - 1),
      js(nghost),
      je(nghost + nj - 1),
      ks(nghost),
      ke(nghost + nk - 1),
      gid(mygid),
      lid(mylid),
      ncells1(ni > 1 ? ni + 2 * nghost + 1 : 2),
      ncells2(nj > 1 ? nj + 2 * nghost + 1 : 2),
      ncells3(nk > 1 ? nk + 2 * nghost + 1 : 2) {
  for (int n = 0; n < 6; ++n) apply_bc[n] = true;
  if (ni == 1) {
    is -= nghost;
    ie -= nghost;
    apply_bc[0] = false;
    apply_bc[1] = false;
  }

  if (nj == 1) {
    js -= nghost;
    je -= nghost;
    apply_bc[2] = false;
    apply_bc[3] = false;
  }

  if (nk == 1) {
    ks -= nghost;
    ke -= nghost;
    apply_bc[4] = false;
    apply_bc[5] = false;
  }
}

BlockIndex::BlockIndex(MeshBlock *pmb)
    : is(pmb->is),
      ie(pmb->ie),
      js(pmb->js),
      je(pmb->je),
      ks(pmb->ks),
      ke(pmb->ke),
      gid(pmb->gid),
      lid(pmb->lid),
      ncells1(pmb->ncells1),
      ncells2(pmb->ncells2),
      ncells3(pmb->ncells3) {
  for (int n = 0; n < 6; ++n) apply_bc[n] = false;
}
