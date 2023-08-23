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

  mass_flux[0].NewAthenaArray(ncells3, ncells2, ncells1 + 1);
  mass_flux[1].NewAthenaArray(ncells3, ncells2 + 1, ncells1);
  mass_flux[2].NewAthenaArray(ncells3 + 1, ncells2, ncells1);

  vsedf[0].NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1 + 1);
  vsedf[1].NewAthenaArray(NCLOUD, ncells3, ncells2 + 1, ncells1);
  vsedf[2].NewAthenaArray(NCLOUD, ncells3 + 1, ncells2, ncells1);

  // hydrodynamic variables
  vsed_.NewAthenaArray(NCLOUD, ncells3, ncells2, ncells1);
  vsed_.ZeroClear();

  hydro_.NewAthenaArray(NCLOUD_HYDRO, NCLOUD, ncells3, ncells2, ncells1);
  hydro_.ZeroClear();

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

  // whether to do sedimentation in X2 or X3 direction
  do_sedimentation_x2_ =
      pin->GetOrAddBoolean("microphysics", "do_sedimentation_x2", false);
  do_sedimentation_x3_ =
      pin->GetOrAddBoolean("microphysics", "do_sedimentation_x3", false);
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
      system->AssembleReactionMatrix(air, time);
      system->EvolveOneStep(&air, time, dt);
    }
}

void Microphysics::UpdateSedimentationVelocity() {
  MeshBlock *pmb = pmy_block_;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  // X1DIR
  for (auto &system : systems_) system->SetSedimentationVelocityX1(vsed_, pmb);

  // interpolation to cell interface
  for (int n = 0; n < NCLOUD; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        for (int i = is + 1; i <= ie; ++i) {
          vsedf[X1DIR](n, k, j, i) =
              interp_cp4(vsed_(n, k, j, i - 2), vsed_(n, k, j, i - 1),
                         vsed_(n, k, j, i), vsed_(n, k, j, i + 1));
        }

        // no sedimentation velocity at the boundary
        vsedf[X1DIR](n, k, j, is) = 0.;
        vsedf[X1DIR](n, k, j, ie + 1) = 0.;
      }

  // X2DIR
  if (pmy_block_->pmy_mesh->f2 && do_sedimentation_x2_) {
    for (auto &system : systems_)
      system->SetSedimentationVelocityX2(vsed_, pmb);

    // interpolation to cell interface
    for (int n = 0; n < NCLOUD; ++n)
      for (int k = ks; k <= ke; ++k) {
        for (int j = js + 1; j <= je; ++j)
          for (int i = is; i <= ie; ++i) {
            vsedf[X2DIR](n, k, j, i) =
                interp_cp4(vsed_(n, k, j - 2, i), vsed_(n, k, j - 1, i),
                           vsed_(n, k, j, i), vsed_(n, k, j + 1, i));
          }

        // no sedimentation velocity at the boundary
        for (int i = is; i <= ie; ++i) {
          vsedf[X2DIR](n, k, js, i) = 0.;
          vsedf[X2DIR](n, k, je + 1, i) = 0.;
        }
      }
  }

  // X3DIR
  if (pmy_block_->pmy_mesh->f3 && do_sedimentation_x3_) {
    for (auto &system : systems_)
      system->SetSedimentationVelocityX3(vsed_, pmb);

    // interpolation to cell interface
    for (int n = 0; n < NCLOUD; ++n) {
      for (int k = ks + 1; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) {
            vsedf[X3DIR](n, k, j, i) =
                interp_cp4(vsed_(n, k - 2, j, i), vsed_(n, k - 1, j, i),
                           vsed_(n, k, j, i), vsed_(n, k + 1, j, i));
          }

      // no sedimentation velocity at the boundary
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          vsedf[X3DIR](n, ks, j, i) = 0.;
          vsedf[X3DIR](n, ke + 1, j, i) = 0.;
        }
    }
  }
}
