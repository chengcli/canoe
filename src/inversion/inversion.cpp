// C/C++ headers
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

// Eigen
#include <Eigen/Core>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// inversion
#include "inversion.hpp"
#include "profile_inversion.hpp"

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name)
    : NamedGroup(name),
      mcmc_initialized_(false),
      init_pos_(nullptr),
      pmy_block_(pmb) {
  Application::Logger app("inversion");
  app->Log("Initialize Inversion");

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

  if (!mcmc_initialized_) {
    app->Error("mcmc chain uninitialized");
  }
  mcmc_save_fits(fname.c_str(), &opts_, &recs_);
}

void Inversion::ResetChain() {
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

AllInversions InversionsFactory::Create(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("inversion");
  app->Log("Create inversion queue");

  std::string str = pin->GetOrAddString("inversion", "tasks", "");
  std::vector<std::string> task_names =
      Vectorize<std::string>(str.c_str(), " ,");

  AllInversions all_fits;
  InversionPtr pfit;

  for (auto p : task_names) {
    if (p == "VLAProfileInversion") {
      pfit = std::make_shared<VLAProfileInversion>(pmb, pin);
    } else if (p == "JunoProfileInversion") {
      pfit = std::make_shared<JunoProfileInversion>(pmb, pin);
    } else if (p == "VLACompositionInversion") {
    } else if (p == "JunoCompositionInversion") {
    } else {
      throw NotFoundError("CreateAllInversions", p);
    }
    all_fits.push_back(std::move(pfit));
  }

  int jl = pmb->js;
  for (auto &q : all_fits) {
    q->InitializePositions();
    q->setX2Indices(jl);
    jl += q->getX2Span();
  }

  app->Log("Number of inversions = " + std::to_string(all_fits.size()));

  return all_fits;
}

namespace InversionHelper {
#define MAX_LINE 512
void read_observation_file(Eigen::VectorXd *target, Eigen::MatrixXd *icov,
                           std::string fname) {
  std::stringstream msg;
  FILE *fp = fopen(fname.c_str(), "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in ProfileInversion::ReadObseravtionFile"
        << std::endl
        << fname << " cannot be opened.";
    ATHENA_ERROR(msg);
  }
  char line[MAX_LINE], *pl;

  int rows;
  // header
  pl = NextLine(line, MAX_LINE, fp);

  // target values
  sscanf(pl, "%d", &rows);
  target->resize(rows);
  icov->resize(rows, rows);

  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    sscanf(pl, "%lf", &(*target)(i));
  }

  // inverse covariance matrix
  char *saveptr;
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok_r(pl, " ", &saveptr);
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &(*icov)(i, j));
      p = strtok_r(NULL, " ", &saveptr);
    }
  }

  fclose(fp);
}
#undef MAX_LINE

void gather_probability(std::vector<Inversion *> const &fitq) {
  auto qlast = fitq[fitq.size() - 1];

  // replace the log probability by the last one
  for (auto q : fitq) {
    int nwalker = q->GetWalkers();

    for (int k = 0; k < nwalker; ++k) {
      q->SetLogProbability(k, qlast->GetLogProbability(k));
    }
  }
}
}  // namespace InversionHelper
