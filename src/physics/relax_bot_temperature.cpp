#include "../communicator/communicator.hpp"
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "physics.hpp"

TaskStatus Physics::RelaxBotTemperature(AthenaArray<Real> &du,
                                        AthenaArray<Real> const &w, Real time,
                                        Real dt) {
  MeshBlock *pmb = pmy_block;
  NeighborBlock const *pbot = pmb->pcomm->findBotNeighbor();
  if (pbot != nullptr) return TaskStatus::success;

  pmb->pdebug->Call("Physics::RelaxBotTemperature");

  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is;
  int js = pmb->js;
  int ks = pmb->ks;
  int ie = pmb->ie;
  int je = pmb->je;
  int ke = pmb->ke;

  Real Rd = pthermo->GetRd();
  Real gamma = pmb->peos->GetGamma();
  Real cv, tem;

  if (Tbot_ < 0.) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        cv = pthermo->getSpecificCv(w.at(k, j, ie));
        tem = pthermo->GetTemp(w.at(k, j, js));
        du(IEN, k, j, is) += dt / tau_Tbot_ * (tem_bot_(k, j) - tem) * cv;
      }
  } else {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        cv = pthermo->getSpecificCv(w.at(k, j, ie));
        tem = pthermo->GetTemp(w.at(k, j, js));
        du(IEN, k, j, is) += dt / tau_Tbot_ * (Tbot_ - tem) * cv;
      }
  }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, is, ie, js, je, ks, ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}
