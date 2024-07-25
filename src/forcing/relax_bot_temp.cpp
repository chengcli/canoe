// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// forcing
#include "forcing.hpp"

RelaxBotTemp::RelaxBotTemp(MeshBlock *pmb, ParameterInput *pin)
    : BotForcing(pmb, 1) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: RelaxBotTemp");

  SetPar("btem", pin->GetReal("forcing", "relax_bot_temp.tem"));
  SetPar("btau", pin->GetReal("forcing", "relax_bot_temp.tau"));
}

void RelaxBotTemp::Initialize(MeshBlock *pmb) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = pmb->phydro->w;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real tem = pthermo->GetTemp(w.at(k, j, pmb->is));
      bot_data_(k, j) = tem;
    }
}

void RelaxBotTemp::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                         Real dt) {
  auto const &w = pmb->phydro->w;
  auto pthermo = Thermodynamics::GetInstance();

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real btem = GetPar<Real>("btem");
  Real btau = GetPar<Real>("btau");

  if (btem < 0.) {  // negative means use the initial condition
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real cv = pthermo->GetCv(w.at(k, j, is));
        Real tem = pthermo->GetTemp(w.at(k, j, is));
        Real rho = w(IDN, k, j, is);
        du(IEN, k, j, is) += dt / btau * (bot_data_(k, j) - tem) * cv * rho;
      }
  } else {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real cv = pthermo->GetCv(w.at(k, j, is));
        Real tem = pthermo->GetTemp(w.at(k, j, is));
        Real rho = w(IDN, k, j, is);
        du(IEN, k, j, is) += dt / btau * (btem - tem) * cv * rho;
      }
  }
}
