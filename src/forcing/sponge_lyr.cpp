// external
#include <application/application.hpp>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/core.h>  // sqr

// forcing
#include "forcing.hpp"

TopSpongeLyr::TopSpongeLyr(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: TopSpongeLyr");

  SetPar("width", pin->GetReal("forcing", "top_sponge_lyr.width"));
  SetPar("tau", pin->GetReal("forcing", "top_sponge_lyr.tau"));
}

BotSpongeLyr::BotSpongeLyr(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: BotSpongeLyr");

  SetPar("width", pin->GetReal("forcing", "bot_sponge_lyr.width"));
  SetPar("tau", pin->GetReal("forcing", "bot_sponge_lyr.tau"));
}

LeftSpongeLyr::LeftSpongeLyr(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: LeftSpongeLyr");

  SetPar("width", pin->GetReal("forcing", "left_sponge_lyr.width"));
  SetPar("tau", pin->GetReal("forcing", "left_sponge_lyr.tau"));
}

RightSpongeLyr::RightSpongeLyr(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: TopSpongeLyr");

  SetPar("width", pin->GetReal("forcing", "right_sponge_lyr.width"));
  SetPar("tau", pin->GetReal("forcing", "right_sponge_lyr.tau"));
}

void TopSpongeLyr::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                         Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  Real x1max = pmb->block_size.x1max;
  Real width = GetPar<Real>("width");
  Real tau = GetPar<Real>("tau");

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->ie; i >= pmb->is; --i) {
        Real eta = (width - (x1max - pcoord->x1f(i + 1))) / width;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau * scale * dt;
      }
}

void BotSpongeLyr::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                         Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  Real x1min = pmb->block_size.x1min;
  Real width = GetPar<Real>("width");
  Real tau = GetPar<Real>("tau");

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real eta = (width - (pcoord->x1f(i) - x1min)) / width;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau * scale * dt;
      }
}

void LeftSpongeLyr::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                          Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  Real x2min = pmb->block_size.x2min;
  Real width = GetPar<Real>("width");
  Real tau = GetPar<Real>("tau");

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real eta = (width - (pcoord->x2f(j) - x2min)) / width;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau * scale * dt;
      }
}

void RightSpongeLyr::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                           Real dt) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  Real x2max = pmb->block_size.x2max;
  Real width = GetPar<Real>("width");
  Real tau = GetPar<Real>("tau");

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->je; j >= pmb->js; --j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real eta = (width - (x2max - pcoord->x2f(j + 1))) / width;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        du(IVX, k, j, i) -=
            w(IDN, k, j, i) * w(IVX, k, j, i) / tau * scale * dt;
        du(IVY, k, j, i) -=
            w(IDN, k, j, i) * w(IVY, k, j, i) / tau * scale * dt;
        du(IVZ, k, j, i) -=
            w(IDN, k, j, i) * w(IVZ, k, j, i) / tau * scale * dt;
      }
}
