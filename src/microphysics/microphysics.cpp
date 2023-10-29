// C/C++
#include <fstream>

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
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>

// microphysics
#include "microphysical_schemes.hpp"
#include "microphysics.hpp"

const std::string Microphysics::input_key = "microphysics_config";

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  if (NCLOUD == 0) return;

  Application::Logger app("microphysics");
  app->Log("Initialize Microphysics");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, 0, NCLOUD);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, 0, NCLOUD);

  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // storage for sedimentation velocity at cell boundary
  vsedf[0].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1 + 1);
  vsedf[1].NewAthenaArray(NCLOUD, ncells3, ncells2 + 1, ncells1);
  vsedf[2].NewAthenaArray(NCLOUD, ncells3 + 1, ncells2, ncells1);

  // storage for cloud mass flux at cell boundary
  vsedf[0].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1 + 1);
  mass_flux[0].NewAthenaArray(ncells3, ncells2, ncells1 + 1);
  mass_flux[1].NewAthenaArray(ncells3, ncells2 + 1, ncells1);
  mass_flux[2].NewAthenaArray(ncells3 + 1, ncells2, ncells1);

  // internal storage for sedimentation velocity at cell center
  vsed_[0].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsed_[0].ZeroClear();
  vsed_[1].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsed_[1].ZeroClear();
  vsed_[2].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsed_[2].ZeroClear();

  // hydro_.NewAthenaArray(NCLOUD_HYDRO, NCLOUD, ncells3, ncells2, ncells1);
  // hydro_.ZeroClear();

  systems_ = MicrophysicalSchemesFactory::Create(pmb, pin);
}

Microphysics::~Microphysics() {
  if (NCLOUD == 0) return;

  Application::Logger app("microphysics");
  app->Log("Destroy Microphysics");
}

// void Microphysics::AddFrictionalHeating(
//     std::vector<AirParcel> &air_column) const {}

void Microphysics::EvolveSystems(AirColumn &air_column, Real time, Real dt) {
  for (auto &system : systems_)
    for (auto &air : air_column) {
      system->AssembleReactionMatrix(air, time);
      system->EvolveOneStep(&air, time, dt);
    }
}

void Microphysics::SetSedimentationVelocityFromConserved(Hydro const *phydro) {
  auto pmb = pmy_block_;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  for (auto &system : systems_) {
    system->SetSedimentationVelocityFromConserved(vsed_, phydro, ks, ke, js, je,
                                                  is, ie);
  }

  // interpolation to cell interface
  for (int n = 0; n < NCLOUD; ++n)
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
  for (int n = 0; n < NCLOUD; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        // no sedimentation velocity at the boundary
        vsedf[X1DIR](n, k, j, is) = 0.;
        vsedf[X1DIR](n, k, j, ie + 1) = 0.;
      }
}

namespace AllTasks {

// hydro tasks should be move into hydro in the future
bool hydro_implicit_correction(MeshBlock *pmb, IntegrationStage stage) {
  return true;
}
bool hydro_calculate_flux(MeshBlock *pmb, IntegrationStage stage) {
  return true;
}

bool microphysics_set_sedimentaton_velocity(MeshBlock *pmb,
                                            IntegrationStage stage) {
  auto scheduler = pmb->pimpl->scheduler;
  if (!scheduler->CheckDone({hydro_implicit_correction})) return false;
  return true;
}

bool microphysics_set_mass_flux(MeshBlock *pmb, IntegrationStage stage) {
  auto scheduler = pmb->pimpl->scheduler;
  if (!scheduler->CheckDone({hydro_calculate_flux})) return false;
  // Riemann Solver in hydro sets the cloud mass flux
  // No need to do anything there
  return true;
}

bool microphysics_evolve_system(MeshBlock *pmb, IntegrationStage stage) {
  auto scheduler = pmb->pimpl->scheduler;
  if (!scheduler->CheckDone({hydro_implicit_correction})) return false;
  return true;
}

bool scalar_claculate_flux(MeshBlock *pmb, IntegrationStage stage) {
  auto scheduler = pmb->pimpl->scheduler;
  if (!scheduler->CheckDone({microphysics_set_mass_flux})) return false;
  return true;
}

}  // namespace AllTasks
