// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include <athena/thermodynamics/thermodynamics.hpp>

// forcing
#include "forcing.hpp"

void RelaxBotComp::RelaxBotComp(MeshBlock *pmb, ParameterInput *pin)
    : HydroForcing(pmb, NVAPOR) {
  Application::Logger app("forcing");
  app->Log("Initialize Forcing: RelaxBotComp");

  for (int i = 1; i <= NVAPOR; ++i) {
    char parname[80];
    char inpname[80];

    snprintf(parname, 80, "bcom%d", i);
    snprintf(inpname, 80, "relax_bot_comp.q%d", i);
    SetPar(parname, pin->GetReal("forcing", inpname));

    snprintf(parname, 80, "btau%d", i);
    snprintf(inpname, 80, "relax_bot_comp.tau%d", i);
    SetPar(parname, pin->GetReal("forcing", inpname));
  }
}

void RelaxBotComp::Initialize(MeshBlock *pmb) {
  auto const &w = pmb->phydro->w;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int n = 0; n < NVAPOR; ++n)
        bot_data_(n, k, j) = w(1 + n, k, j, pmb->is);
}

void RelaxBotComp::Apply(AthenaArray<Real> &du, MeshBlock *pmb, Real time,
                         Real dt) {
  auto pthermo = Thermodynamics::GetInstance();
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real gammad = pthermo->GetGammadRef();
  Real cvd = pthermo->GetRd() / (gammad - 1.);

  for (int n = 0; n < NVAPOR; ++n) {
    Real bcom = GetPar<Real>("bcom" + std::to_string(n));
    Real cv_ratio = pthermo->GetCvRatioMass(n);

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real tem = pthermo->GetTemp(pmb, k, j, is);
        Real v1 = w(IVX, k, j, is);
        Real v2 = w(IVY, k, j, is);
        Real v3 = w(IVZ, k, j, is);
        Real cv0 = pthermo

            if (bcom < 0.) {  // negative means use the initial condition
          Real dq = bot_data_(n, k, j) - w(1 + n, k, j, is);
        }
        else {
          Real dq = bcom - w(1 + n, k, j, is);
        }

        du(IDN, k, j, is) -= dt / btau * dq * rho;
        du(1 + n, k, j, is) += dt / btau * dq * rho;

        du(IVX, k, j, is) += dt / btau * v1 * rho * dq du(IVY, k, j, is) +=
            dt / btau * v2 * rho * dq du(IVZ, k, j, is) +=
            dt / btau * v3 * rho *
            dq

                du(IEN, k, j, is) +=
            dt / btau * dq * (1. - cv_ratio) * rho * tem * cvd;
      }
  }
}
