//! \brief Implementation of the Microphysics class
//! \file microphysics.cpp

// C/C++
#include <fstream>

// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/reconstruct/interpolation.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// microphysics
#include "microphysics.hpp"

const std::string Microphysics::input_key = "microphysics_config";

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  if (NMASS == 0) return;

  std::string fname = pin->GetString("problem", input_key);

  Application::Logger app("microphysics");
  app->Log("Initialize Microphysics");

  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // storage for sedimentation velocity at cell boundary
  vsedf[0].NewAthenaArray(NMASS, ncells3, ncells2, ncells1 + 1);
  vsedf[1].NewAthenaArray(NMASS, ncells3, ncells2 + 1, ncells1);
  vsedf[2].NewAthenaArray(NMASS, ncells3 + 1, ncells2, ncells1);

  // internal storage for sedimentation velocity at cell center
  vsed_[0].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[0].ZeroClear();
  vsed_[1].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[1].ZeroClear();
  vsed_[2].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[2].ZeroClear();

  systems_ = MicrophysicalSchemesFactory::Create(pmb, pin);
}

Microphysics::~Microphysics() {
  if (NMASS == 0) return;

  Application::Logger app("microphysics");
  app->Log("Destroy Microphysics");
}

// void Microphysics::AddFrictionalHeating(
//     std::vector<AirParcel> &air_column) const {}

void Microphysics::EvolveSystems(T u, T s, Real time, Real dt) {
  for (auto &sys : systems_) {
    sys->EvolveOneStep(u, s, time, dt);
  }
}

void Microphysics::SetVsedFromConserved(Hydro const *phydro) {
  auto pmb = pmy_block_;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  for (auto &sys : systems_) {
    sys->SetVsedFromConserved(vsed_, phydro, ks, ke, js, je, is, ie);
  }

  // interpolation to cell interface
  for (int n = 0; n < NMASS; ++n)
    for (int k = ks; k <= ke + 1; ++k)
      for (int j = js; j <= je + 1; ++j)
        for (int i = is; i <= ie + 1; ++i) {
          vsedf[X1DIR](n, k, j, i) = interp_cp4(
              vsed_[X1DIR](n, k, j, i - 2), vsed_[X1DIR](n, k, j, i - 1),
              vsed_[X1DIR](n, k, j, i), vsed_[X1DIR](n, k, j, i + 1));

          vsedf[X2DIR](n, k, j, i) = interp_cp4(
              vsed_[X2DIR](n, k, j - 2, i), vsed_[X2DIR](n, k, j - 1, i),
              vsed_[X2DIR](n, k, j, i), vsed_[X2DIR](n, k, j + 1, i));

          vsedf[X3DIR](n, k, j, i) = interp_cp4(
              vsed_[X3DIR](n, k - 2, j, i), vsed_[X3DIR](n, k - 1, j, i),
              vsed_[X3DIR](n, k, j, i), vsed_[X3DIR](n, k + 1, j, i));
        }

  // fix boundary condition (TODO)
  for (int n = 0; n < NMASS; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        // no sedimentation velocity at the boundary
        vsedf[X1DIR](n, k, j, is) = 0.;
        vsedf[X1DIR](n, k, j, ie + 1) = 0.;
      }
}
