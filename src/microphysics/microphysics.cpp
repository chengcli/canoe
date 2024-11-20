//! \brief Implementation of the Microphysics class
//! \file microphysics.cpp

// C/C++
#include <fstream>

// athena
#include <athena/athena.hpp>
#include <athena/hydro/hydro.hpp>
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
#include <snap/athena_torch.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// utils
#include <utils/vectorize.hpp>

// microphysics
#include "microphysical_schemes.hpp"
#include "microphysics.hpp"

const std::string Microphysics::input_key = "microphysics_config";

enum { NMASS = NCLOUD + NPRECIP };

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb) {
  if (NMASS == 0) return;

  Application::Logger app("microphysics");
  app->Log("Initialize Microphysics");

  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // storage for sedimentation velocity at cell boundary
  vsedf[0].NewAthenaArray(NMASS, ncells3, ncells2, ncells1 + 1);
  vsedf[1].NewAthenaArray(NMASS, ncells3, ncells2 + 1, ncells1);
  vsedf[2].NewAthenaArray(NMASS, ncells3 + 1, ncells2, ncells1);

  // storage for cloud mass flux at cell boundary
  mass_flux[0].NewAthenaArray(NMASS, ncells3, ncells2, ncells1 + 1);
  mass_flux[1].NewAthenaArray(NMASS, ncells3, ncells2 + 1, ncells1);
  mass_flux[2].NewAthenaArray(NMASS, ncells3 + 1, ncells2, ncells1);

  // internal storage for sedimentation velocity at cell center
  vsed_[0].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[0].ZeroClear();
  vsed_[1].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[1].ZeroClear();
  vsed_[2].NewAthenaArray(NMASS, ncells3, ncells2, ncells1);
  vsed_[2].ZeroClear();

  systems_ = MicrophysicalSchemesFactory::Create(pmb, pin);

  // set up sedimentation options
  auto str = pin->GetOrAddString("microphysics", "particle_radius", "0.");
  sed_opts_.radius() = Vectorize<double>(str.c_str(), " ,");

  str = pin->GetOrAddString("microphysics", "particle_density", "0.");
  sed_opts_.density() = Vectorize<double>(str.c_str(), " ,");

  auto vsed1 = pin->GetOrAddReal("microphysics", "particle_vsed1", 0.);
  auto vsed2 = pin->GetOrAddReal("microphysics", "particle_vsed2", 0.);

  if (vsed1 != 0.0) {
    sed_opts_.const_vsed().push_back(vsed1);
    if (vsed2 != 0.0) sed_opts_.const_vsed().push_back(vsed2);
  } else if (vsed2 != 0.0) {
    sed_opts_.const_vsed().push_back(0.);
    sed_opts_.const_vsed().push_back(vsed2);
  }

  sed_opts_.gravity() = pin->GetOrAddReal("hydro", "grav_acc1", 0.0);
}

Microphysics::~Microphysics() {
  if (NMASS == 0) return;

  Application::Logger app("microphysics");
  app->Log("Destroy Microphysics");
}

// void Microphysics::AddFrictionalHeating(
//     std::vector<AirParcel> &air_column) const {}

void Microphysics::EvolveSystems(AirColumn &ac, Real time, Real dt) {
  for (auto &system : systems_)
    for (auto &air : ac) {
      system->AssembleReactionMatrix(air, time);
      system->EvolveOneStep(&air, time, dt);
    }
}

void Microphysics::SetVsedFromConserved(Hydro const *phydro) {
  if (NMASS == 0) return;

  auto pmb = pmy_block_;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  // for (auto &system : systems_) {
  //   system->SetVsedFromConserved(vsed_, phydro, ks, ke, js, je, is, ie);
  // }

  // Following this article
  // https://stackoverflow.com/questions/59676983/in-torch-c-api-how-to-write-to-the-internal-data-of-a-tensor-fastly

  // set temperature
  auto pthermo = Thermodynamics::GetInstance();
  torch::Tensor temp =
      torch::zeros({pmb->ncells3, pmb->ncells2, pmb->ncells1}, torch::kDouble);
  auto accessor = temp.accessor<double, 3>();

  for (int k = 0; k < pmb->ncells3; ++k)
    for (int j = 0; j < pmb->ncells2; ++j)
      for (int i = 0; i < pmb->ncells1; ++i) {
        accessor[k][j][i] = pthermo->GetTemp(phydro->w.at(k, j, i));
      }

  auto sed = Sedimentation(sed_opts_);
  sed->to(torch::kCPU, torch::kDouble);

  auto diag = std::make_shared<SharedData::element_type>();
  (*diag)["temperature"] =
      std::async(std::launch::async, [&]() { return temp; }).share();
  sed->set_shared_data(diag);

  // calculate sedimentation velocity
  auto vel = sed->forward(to_torch(phydro->w));

  vsed_[X1DIR] = to_athena(vel);

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
