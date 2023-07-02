// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"
#include "inversion_helper.hpp"
#include "profile_inversion.hpp"

InversionQueue Inversion::NewInversionQueue(MeshBlock *pmb,
                                            ParameterInput *pin) {
  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names =
      Vectorize<std::string>(str.c_str(), " ,");

  InversionQueue fitq;
  Application::Logger app("inversion");
  app->Log("Create inversion queue");

  InversionPtr pfit;

  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = std::make_unique<VLAProfileInversion>(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = std::make_unique<JunoProfileInversion>(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      throw NotFoundError("NewInversionQueue", p);
    }
    fitq.push_back(std::move(pfit));
  }

  int jl = pmb->js;
  for (auto &q : fitq) {
    q->InitializePositions();
    q->setX2Indices(jl);
    jl += q->getX2Span();
  }

  app->Log("Number of inversions = " + std::to_string(fitq.size()));

  return fitq;
}

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name)
    : name_(name),
      mcmc_initialized_(false),
      init_pos_(nullptr),
      pmy_block_(pmb) {
  Application::Logger app("inversion");
  app->Log("Initialize Inversion");
  std::stringstream msg;

  opts_.a = pin->GetOrAddReal("inversion", "stretch", 2.);
  opts_.p = pin->GetOrAddInteger("inversion", "walk", 4);
  opts_.print = pin->GetOrAddInteger("inversion", "print", 100);

#ifdef MPI_PARALLEL
  opts_.mpi_comm = MPI_COMM_WORLD;
#endif

  snprintf(opts_.logfile, sizeof(opts_.logfile), "%s",
           pin->GetOrAddString("inversion", name + ".logfile",
                               name + "_inversion.log")
               .c_str());

  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "none");
  if (obsfile != "none") {
    read_observation_file(&target_, &icov_, obsfile.c_str());
    // app->Log("target = ", target_.transpose());
    // app->Log("inverse covariance matrx = ", icov_);
  } else {
    target_.resize(1);
    icov_.resize(1, 1);
  }

  // fit differential
  fit_differential_ = pin->GetOrAddBoolean("inversion", "differential", false);
  // app->Log("fit differential = ", fit_differential_);
}

Inversion::~Inversion() {
  if (mcmc_initialized_) {
    mcmc_free(&recs_);
    delete[] zz_;
    delete[] par_;
  }

  if (init_pos_ != nullptr) FreeCArray(init_pos_);
}

void Inversion::InitializeChain(int nstep, int nwalker, int ndim, int nvalue) {
  // mcmc_alloc(&recs_, pmy_block->pmy_mesh->nlim+1, nwalker, ndim, nvalue);
  mcmc_alloc(&recs_, nstep, nwalker, ndim, nvalue);
  mcmc_initialized_ = true;
  zz_ = new Real[nwalker];
  par_ = new Real[ndim];
}

void Inversion::MakeMCMCOutputs(std::string fname) {
  Application::Logger app("inversion");
  app->Log("Make MCMC Outputs");

  std::stringstream msg;
  if (!mcmc_initialized_) {
    app->Error("mcmc chain uninitialized");
  }
  mcmc_save_fits(fname.c_str(), &opts_, &recs_);
}

void Inversion::ResetChain() {
  std::stringstream msg;
  Application::Logger app("inversion");
  app->Log("Reset MCMC Chain");

  if (!mcmc_initialized_) {
    app->Error("mcmc chain uninitialized");
  }

  int cur = recs_.cur;
  if (cur == 0) return;

  // copy the last state into the first state
  for (int k = 0; k < recs_.nwalker; ++k) {
    for (int d = 0; d < recs_.ndim; ++d)
      recs_.par[0][k][d] = recs_.par[cur - 1][k][d];
    for (int d = 0; d < recs_.nvalue; ++d)
      recs_.val[0][k][d] = recs_.val[cur - 1][k][d];
    recs_.lnp[0][k] = recs_.lnp[cur - 1][k];
    recs_.newstate[0][k] = recs_.newstate[cur - 1][k];
  }

  recs_.reset += cur - 1;
  recs_.cur = 1;

  int nstep = recs_.nstep - 1;
  int nwalker = recs_.nwalker;
  int ndim = recs_.ndim;
  int nvalue = recs_.nvalue;

  memset(recs_.par[1][0], 0, nstep * nwalker * ndim * sizeof(double));
  memset(recs_.val[1][0], 0, nstep * nwalker * nvalue * sizeof(double));
  memset(recs_.lnp[1], 0, nstep * nwalker * sizeof(double));
  memset(recs_.newstate[1], 0, nstep * nwalker * sizeof(int));
}
