// Athena++ headers
#include "../communicator/communicator.hpp"
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "physics.hpp"

TaskStatus Physics::TopCooling(AthenaArray<Real> &du,
                               AthenaArray<Real> const &w, Real time, Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *ptop = pmb->pcomm->findTopNeighbor();
  if (ptop != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::TopCooling");

  Coordinates *pcoord = pmb->pcoord;
  int is = pmb->is;
  int js = pmb->js;
  int ks = pmb->ks;
  int ie = pmb->ie;
  int je = pmb->je;
  int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      // Real cv = pthermo->GetMeanCv(w.at(k,j,ie));
      // du(IEN,k,j,ie) -= dt*dTdt_*w(IDN,k,j,ie)*cv;
      du(IEN, k, j, ie) += dt * flux_top_ * pcoord->GetFace1Area(k, j, ie + 1) /
                           pcoord->GetCellVolume(k, j, ie);
    }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, is, ie, js, je, ks, ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

TaskStatus Physics::BotHeating(AthenaArray<Real> &du,
                               AthenaArray<Real> const &w, Real time, Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *pbot = pmb->pcomm->findBotNeighbor();
  if (pbot != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::BotHeating");

  Coordinates *pcoord = pmb->pcoord;
  int is = pmb->is;
  int js = pmb->js;
  int ks = pmb->ks;
  int ie = pmb->ie;
  int je = pmb->je;
  int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      du(IEN, k, j, is) += dt * flux_bot_ * pcoord->GetFace1Area(k, j, is) /
                           pcoord->GetCellVolume(k, j, is);
    }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, is, ie, js, je, ks, ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}
