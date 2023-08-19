// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// microphysics
#include "microphysics.hpp"

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  if (NCLOUD == 0) return;

  Application::Logger app("microphysics");
  app->Log("Initialize Microphysics");

  w.InitWithShallowSlice(pmb->pscalars->r, 4, 0, NCLOUD);
  u.InitWithShallowSlice(pmb->pscalars->s, 4, 0, NCLOUD);

  // hydrodynamic variables
  vsed_.NewAthenaArray(NCLOUD, pmb->ncells1);
  vsedf_.NewAthenaArray(NCLOUD, pmb->ncells1 + 1);
  hydro_.NewAthenaArray(NCLOUD_HYDRO, NCLOUD, pmb->ncells3, pmb->ncells2,
                        pmb->ncells1);
}

Microphysics::~Microphysics() {
  if (NCLOUD == 0) return;

  Application::Logger app("microphysics");
  app->Log("Destroy Microphysics");
}

void Microphysics::AddFrictionalHeating(std::vector<AirParcel> &air_column) {}

void Microphysics::EvolveSystems(std::vector<AirParcel> &air_column, Real time,
                                 Real dt) {
  for (auto &system : systems_)
    for (auto &air : air_column) {
      system->AssembleReactionMatrix(system->GetRatePtr(),
                                     system->GetJacobianPtr(), air, time);
      system->EvolveOneStep(&air, time, dt);
    }
}

void Microphysics::SetSedimentationVelocity(int k, int j, int il, int iu) {
  for (auto &system : systems_) {
    system->SetSedimentationVelocity(vsed_, k, j, il, iu);
  }

  // interpolation to cell interface
}

void Microphysics::AddSedimentationFlux(AthenaArray<Real> &sflx, int k, int j,
                                        int il, int iu) {}
