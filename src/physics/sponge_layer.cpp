// Athena++ headers
#include "../communicator/communicator.hpp"
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "../math/core.h"  // sqr
#include "../mesh/mesh.hpp"
#include "physics.hpp"

TaskStatus Physics::TopSpongeLayer(AthenaArray<Real> &du,
                                   AthenaArray<Real> const &w, Real time,
                                   Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *ptop = pmb->pcomm->findTopNeighbor();
  if (ptop != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::TopSpongerLayer");
  Coordinates *pcoord = pmb->pcoord;
  Real x1max = pmb->block_size.x1max;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->ie; i >= pmb->is; --i) {
        Real eta = (width_top_ - (x1max - pcoord->x1f(i + 1))) / width_top_;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau_top_ * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau_top_ * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau_top_ * scale * dt;
      }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

TaskStatus Physics::BotSpongeLayer(AthenaArray<Real> &du,
                                   AthenaArray<Real> const &w, Real time,
                                   Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *pbot = pmb->pcomm->findBotNeighbor();
  if (pbot != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::BotSpongeLayer");
  Coordinates *pcoord = pmb->pcoord;
  Real x1min = pmb->block_size.x1min;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real eta = (width_bot_ - (pcoord->x1f(i) - x1min)) / width_bot_;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau_bot_ * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau_bot_ * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau_bot_ * scale * dt;
      }

#if (DEBUG_LEVEL > 1)
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

TaskStatus Physics::LeftSpongeLayer(AthenaArray<Real> &du,
                                    AthenaArray<Real> const &w, Real time,
                                    Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *pleft = pmb->pcomm->findLeftNeighbor();
  if (pleft != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::LeftSpongerLayer");
  Coordinates *pcoord = pmb->pcoord;
  Real x2min = pmb->block_size.x2min;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real eta = (width_left_ - (pcoord->x2f(j) - x2min)) / width_left_;
      if (eta < 0) break;
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau_left_ * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau_left_ * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau_left_ * scale * dt;
      }
    }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

TaskStatus Physics::RightSpongeLayer(AthenaArray<Real> &du,
                                     AthenaArray<Real> const &w, Real time,
                                     Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *pright = pmb->pcomm->findRightNeighbor();
  if (pright != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::RightSpongeLayer");
  Coordinates *pcoord = pmb->pcoord;
  Real x2max = pmb->block_size.x2max;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->je; j <= pmb->js; --j) {
      Real eta = (width_right_ - (x2max - pcoord->x2f(j + 1))) / width_right_;
      if (eta < 0) break;
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau_right_ * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau_right_ * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau_right_ * scale * dt;
      }
    }

#if (DEBUG_LEVEL > 1)
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}
