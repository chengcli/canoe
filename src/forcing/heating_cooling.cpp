// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// forcing
#include "forcingt.hpp"

void TopCooling::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                       Real dt) {
  Coordinates *pcoord = pmb->pcoord;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real flux = GetPar<Real>("flux");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      // Real cv = pthermo->GetMeanCv(w.at(k,j,ie));
      // du(IEN,k,j,ie) -= dt*dTdt_*w(IDN,k,j,ie)*cv;
      du(IEN, k, j, ie) += dt * flux * pcoord->GetFace1Area(k, j, ie + 1) /
                           pcoord->GetCellVolume(k, j, ie);
    }
}

void BotHeating::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                       Real dt) {
  Coordinates *pcoord = pmb->pcoord;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real flux = GetPar<Real>("flux");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      du(IEN, k, j, is) += dt * flux * pcoord->GetFace1Area(k, j, is) /
                           pcoord->GetCellVolume(k, j, is);
    }
}

void BodyHeating::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                        Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;
  auto pthermo = Thermodynamics::GetInstance();

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real flux = GetPar<Real>("dTdt");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real cv = pthermo->GetMeanCv(w.at(k, j, ie));
        du(IEN, k, j, ie) += dt * dTdt * w(IDN, k, j, ie) * cv;
      }
}
