// C/C++
#include <fstream>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/scalars/scalars.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// microphysics
#include "microphysical_scheme.hpp"
#include "microphysics.hpp"

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

  // hydrodynamic variables
  vsed_.NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsedf_.NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1 + 1);
  hydro_.NewAthenaArray(NCLOUD_HYDRO, NCLOUD, ncells3, ncells2, ncells1);

  // load all microphysics systems
  std::string key = "microphysics_config";

  if (pin->DoesParameterExist("chemistry", key)) {
    std::string filename = pin->GetString("chemistry", key);
    std::ifstream stream(filename);
    if (stream.good() == false) {
      app->Error("Cannot open microphysics config file: " + filename);
    }

    YAML::Node node = YAML::Load(stream);
    if (!node["microphysics"]) {
      throw NotFoundError("Microphysics", "microphysics");
    }

    for (auto sys : node["microphysics"]) {
      std::string name = sys.as<std::string>();
      std::string scheme = node[name]["scheme"].as<std::string>();
      if (scheme == "Kessler94") {
        auto p = std::make_unique<Kessler94>(name, node[name]);
        systems_.push_back(std::move(p));
      } else {
        throw NotFoundError("Microphysics", scheme);
      }
    }
  }
}

Microphysics::~Microphysics() {
  if (NCLOUD == 0) return;

  Application::Logger app("microphysics");
  app->Log("Destroy Microphysics");
}

void Microphysics::AddFrictionalHeating(
    std::vector<AirParcel> &air_column) const {}

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
                                        int il, int iu) const {}
