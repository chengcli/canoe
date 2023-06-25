// C/C++
#include <memory>
#include <string>

// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// harp
#include "concentration_inversion.hpp"

ConcentrationInversion::~ConcentrationInversion() {}

ConcentrationInversion::ConcentrationInversion(MeshBlock *pmb,
                                               ParameterInput *pin,
                                               std::string name)
    : Inversion(pmb, pin, name) {
  Application::Logger app("inversion");
  app->Log("Initializing ConcentrationInversion");
  char buf[80];

  // species id
  idx_ =
      Vectorize<int>(pin->GetString("inversion", name + ".variables").c_str());

  // read in prior
  for (auto m : idx_) {
    if (m == IDN) {  // change temperature
      Xstd_[IDN] = pin->GetReal("inversion", name + ".tem.std");
      app->Log(name + "::temperature std = " + std::to_string(Xstd_[IDN]));
    } else {
      Xstd_[m] = pin->GetReal("inversion", name + ".qvapor" +
                                               std::to_string(m) + ".std.gkg") /
                 1.E3;
      snprintf(buf, sizeof(buf), "%s::vapor %d standard deviation",
               name.c_str(), m);
      app->Log(buf + std::to_string(Xstd_[m]));
    }
  }

  int ndim = idx_.size();
  app->Log(name + "::number of input dimension = " + std::to_string(ndim));

  // output dimension
  int nvalue = target_.size();
  app->Log("number of output dimension = " + std::to_string(nvalue));

  // number of walkers
  int nwalker = pmb->block_size.nx3;
  app->Log("walkers per block = " + std::to_string(nwalker));
  app->Log("total number of walkers = " +
           std::to_string(pmb->pmy_mesh->mesh_size.nx3));
  if ((nwalker < 2) && pmb->pmy_mesh->nlim > 0) {
    app->Error("nwalker (nx3) must be at least 2");
  }

  // initialize mcmc chain
  InitializeChain(pmb->pmy_mesh->nlim + 1, nwalker, ndim, nvalue);
}

void ConcentrationInversion::InitializePositions() {
  int nwalker = GetWalkers();
  int ndim = GetDims();
  Application::Logger app("inversion");

  // initialize random positions
  app->Log("Initializing random positions for walkers");

  unsigned int seed = time(NULL) + Globals::my_rank;
  NewCArray(init_pos_, nwalker, ndim);

  for (int p = 0; p < nwalker; ++p) {
    for (size_t n = 0; n < idx_.size(); ++n) {
      int m = idx_[n];
      init_pos_[p][n] = (1. * rand_r(&seed) / RAND_MAX - 0.5) * Xstd_[m];
    }
  }
}

void ConcentrationInversion::UpdateConcentration(Hydro *phydro, Real *Xp, int k,
                                                 int jl, int ju) const {
  Application::Logger app("inversion");
  app->Log("UpdateConcentration");

  // int is = pblock_->is, ie = pblock_->ie;
}
