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

  mass_flux_[0].NewAthenaArray(ncells3, ncells2, ncells1 + 1);
  mass_flux_[1].NewAthenaArray(ncells3, ncells2 + 1, ncells1);
  mass_flux_[2].NewAthenaArray(ncells3 + 1, ncells2, ncells1);

  vsedf_[0].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1 + 1);
  vsedf_[1].NewAthenaArray(NCLOUD, ncells3, ncells2 + 1, ncells1);
  vsedf_[2].NewAthenaArray(NCLOUD, ncells3 + 1, ncells2, ncells1);

  // temporary storage for sedimentation velocity
  vsed_.NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsed_.ZeroClear();

  // hydro_.NewAthenaArray(NCLOUD_HYDRO, NCLOUD, ncells3, ncells2, ncells1);
  // hydro_.ZeroClear();

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

// void Microphysics::AddFrictionalHeating(
//     std::vector<AirParcel> &air_column) const {}

void Microphysics::EvolveSystems(AirColumn &air_column, Real time, Real dt) {
  for (auto &system : systems_)
    for (auto &air : air_column) {
      system->AssembleReactionMatrix(air, time);
      system->EvolveOneStep(&air, time, dt);
    }
}

void Microphysics::UpdateSedimentationVelocityFromConserved() {
  MeshBlock *pmb = pmy_block_;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  // X1DIR
  for (auto &system : systems_)
    system->SetSedimentationVelocityFromConserved(pmb->phydro, ks, ke, js, je,
                                                  is, ie);

  // interpolation to cell interface
  for (int n = 0; n < NCLOUD; ++n)
    for (int k = ks; k <= ke + 1; ++k)
      for (int j = js; j <= je + 1; ++j)
        for (int i = is; i <= ie + 1; ++i) {
          vsedf_[X1DIR](n, k, j, i) = interp_cp4(
              vsed_[X1DIR](n, k, j, i - 2), vsed_[X1DIR](n, k, j, i - 1),
              vsed_[X1DIR](n, k, j, i), vsed_[X1DIR](n, k, j, i + 1));

          vsedf_[X2DIR](n, k, j, i) = interp_cp4(
              vsed_[X2DIR](n, k, j - 2, i), vsed_[X2DIR](n, k, j - 1, i),
              vsed_[X2DIR](n, k, j, i), vsed_[X2DIR](n, k, j + 1, i));

          vsedf_[X3DIR](n, k, j, i) = interp_cp4(
              vsed_[X3DIR](n, k - 2, j, i), vsed_[X3DIR](n, k - 1, j, i),
              vsed_[X3DIR](n, k, j, i), vsed_[X3DIR](n, k + 1, j, i));
        }

  // fix boundary condition (TODO)
  for (int n = 0; n < NCLOUD; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        // no sedimentation velocity at the boundary
        vsedf_[X1DIR](n, k, j, is) = 0.;
        vsedf_[X1DIR](n, k, j, ie + 1) = 0.;
      }
}
}
