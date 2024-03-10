// external
#include <application/application.hpp>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// forcing
#include "forcing.hpp"

RelaxBotVelo::RelaxBotVelo(MeshBlock *pmb, ParameterInput *pin)
    : BotForcing(pmb, 3) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: RelaxBotVelo");

  for (int i = 1; i <= NVAPOR; ++i) {
    char parname[80];
    char inpname[80];

    snprintf(parname, 80, "bv%d", i);
    snprintf(inpname, 80, "relax_bot_velo.v%d", i);
    SetPar(parname, pin->GetReal("forcing", inpname));

    snprintf(parname, 80, "btau%d", i);
    snprintf(inpname, 80, "relax_bot_velo.tau%d", i);
    SetPar(parname, pin->GetReal("forcing", inpname));

    snprintf(parname, 80, "use_v%d_init", i);
    snprintf(inpname, 80, "relax_bot_velo.use_v%d_init", i);
    SetPar(parname, pin->GetOrAddInteger("forcing", inpname, 0));
  }
}

void RelaxBotVelo::Initialize(MeshBlock *pmb) {
  auto const &w = pmb->phydro->w;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      bot_data_(0, k, j) = w(IVX, k, j, pmb->is);
      bot_data_(1, k, j) = w(IVY, k, j, pmb->is);
      bot_data_(2, k, j) = w(IVZ, k, j, pmb->is);
    }
}

void RelaxBotVelo::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                         Real dt) {
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real btau1 = GetPar<Real>("btau1");
  Real btau2 = GetPar<Real>("btau2");
  Real btau3 = GetPar<Real>("btau3");

  Real bv1 = GetPar<Real>("bv1");
  Real bv2 = GetPar<Real>("bv2");
  Real bv3 = GetPar<Real>("bv3");

  int use_v1_init = GetPar<int>("use_v1_init");
  int use_v2_init = GetPar<int>("use_v2_init");
  int use_v3_init = GetPar<int>("use_v3_init");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real rho = w(IDN, k, j, is);
      Real v1 = w(IVX, k, j, is);
      if (use_v1_init) {
        du(IVX, k, j, is) += dt / btau1 * (bot_data_(0, k, j) - v1) * rho;
      } else {
        du(IVX, k, j, is) += dt / btau1 * (bv1 - v1) * rho;
      }

      Real v2 = w(IVY, k, j, is);
      if (use_v2_init) {
        du(IVY, k, j, is) += dt / btau2 * (bot_data_(1, k, j) - v2) * rho;
      } else {
        du(IVY, k, j, is) += dt / btau2 * (bv2 - v2) * rho;
      }

      Real v3 = w(IVZ, k, j, is);
      if (use_v3_init) {
        du(IVZ, k, j, is) += dt / btau3 * (bot_data_(2, k, j) - v3) * rho;
      } else {
        du(IVZ, k, j, is) += dt / btau3 * (bv3 - v3) * rho;
      }
    }
}
