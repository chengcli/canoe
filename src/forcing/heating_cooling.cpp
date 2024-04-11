// external
#include <application/application.hpp>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// forcing
#include "forcing.hpp"

TopCooling::TopCooling(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: TopCooling");

  SetPar("flux", pin->GetReal("forcing", "top_cooling.flux"));
  srand(Globals::my_rank + time(0));
}

BotHeating::BotHeating(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: BotHeating");

  SetPar("flux", pin->GetReal("forcing", "bot_heating.flux"));
  srand(Globals::my_rank + time(0));
}

BodyHeating::BodyHeating(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: BodyHeating");

  SetPar("dTdt", pin->GetReal("forcing", "body_heating.dTdt.Kday") / 86400.);
  SetPar("pmin", pin->GetReal("forcing", "body_heating.pmin"));
  SetPar("pmax", pin->GetReal("forcing", "body_heating.pmax"));
  srand(Globals::my_rank + time(0));
}

void TopCooling::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                       Real dt) {
  Coordinates *pcoord = pmb->pcoord;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real flux = GetPar<Real>("flux");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      du(IEN, k, j, ie) +=
          dt * flux * (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX)) *
          pcoord->GetFace1Area(k, j, ie + 1) / pcoord->GetCellVolume(k, j, ie);
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
      du(IEN, k, j, is) +=
          dt * flux * (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX)) *
          pcoord->GetFace1Area(k, j, is) / pcoord->GetCellVolume(k, j, is);
    }
}

void BodyHeating::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                        Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;
  auto pthermo = Thermodynamics::GetInstance();

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real dTdt = GetPar<Real>("dTdt");
  Real pmin = GetPar<Real>("pmin");
  Real pmax = GetPar<Real>("pmax");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        if (w(IDN, k, j, i) < pmin || w(IDN, k, j, i) > pmax) continue;
        Real cv = pthermo->GetCvMass(pmb, k, j, i);
        du(IEN, k, j, ie) += dt * dTdt * w(IDN, k, j, i) * cv;
      }
}
