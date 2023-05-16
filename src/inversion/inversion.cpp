// C/C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

// Athena++ headers
#include <parameter_input.hpp>
#include <mesh/mesh.hpp>
#include <coordinates/coordinates.hpp>

// harp headers
#include "../utils/ndarrays.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/block_index.hpp"
#include "../mesh/meshblock_impl.hpp"
#include "../utils/sentinelq.hpp"
#include "inversion.hpp"

extern std::unique_ptr<Debugger> pdebug;

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name):
  name_(name), mcmc_initialized_(false), init_pos_(nullptr),
  pcoord_(pmb->pcoord)
{
  pdebug->Enter("Inversion");
  std::stringstream msg;

  opts_.a = pin->GetOrAddReal("inversion", "stretch", 2.);
  opts_.p = pin->GetOrAddInteger("inversion", "walk", 4);
  opts_.print = pin->GetOrAddInteger("inversion", "print", 100);

#ifdef MPI_PARALLEL
  opts_.mpi_comm = MPI_COMM_WORLD;
#endif

  strcpy(opts_.logfile, pin->GetOrAddString("inversion", name + ".logfile",
    name + "_inversion.log").c_str());

  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "none");
  if (obsfile != "none") {
    read_observation_file(target_, icov_, obsfile.c_str());
    pdebug->Message("target", target_.transpose());
    pdebug->Message("inverse covariance matrx", icov_);
  } else {
    target_.resize(1);
    icov_.resize(1,1);
  }

  // fit differential
  fit_differential_ = pin->GetOrAddBoolean("inversion", "differential", false);
  pdebug->Message("fit differential", fit_differential_);

  pdebug->Leave();
}

Inversion::~Inversion()
{
  if (mcmc_initialized_) {
    mcmc_free(&recs_);
    delete[] zz_;
    delete[] par_;
  }

  if (init_pos_ != nullptr)
    FreeCArray(init_pos_);
}

void Inversion::InitializeChain(int nstep, int nwalker, int ndim, int nvalue)
{
	//mcmc_alloc(&recs_, pmy_block->pmy_mesh->nlim+1, nwalker, ndim, nvalue);
	mcmc_alloc(&recs_, nstep, nwalker, ndim, nvalue);
	mcmc_initialized_ = true;
  zz_ = new Real [nwalker];
  par_ = new Real [ndim];
}

void Inversion::MakeMCMCOutputs(std::string fname)
{
  std::stringstream msg;
  if (!mcmc_initialized_) {
    Debugger::Fatal("MakeMCMCOutputs", "mcmc chain uninitialized");
  }
  mcmc_save_fits(fname.c_str(), &opts_, &recs_);
}

void Inversion::ResetChain()
{
  std::stringstream msg;
  if (!mcmc_initialized_) {
    Debugger::Fatal("ResetChain", "mcmc chain uninitialized");
  }

  int cur = recs_.cur;
  if (cur == 0) return;

  // copy the last state into the first state
  for (int k = 0; k < recs_.nwalker; ++k) {
    for (int d = 0; d < recs_.ndim; ++d)
      recs_.par[0][k][d] = recs_.par[cur-1][k][d];
    for (int d = 0; d < recs_.nvalue; ++d)
      recs_.val[0][k][d] = recs_.val[cur-1][k][d];
    recs_.lnp[0][k] = recs_.lnp[cur-1][k];
    recs_.newstate[0][k] = recs_.newstate[cur-1][k];
  }

  recs_.reset += cur - 1;
  recs_.cur = 1;

  int nstep = recs_.nstep - 1;
  int nwalker = recs_.nwalker;
  int ndim = recs_.ndim;
  int nvalue = recs_.nvalue;

  memset(recs_.par[1][0], 0, nstep*nwalker*ndim*sizeof(double));
  memset(recs_.val[1][0], 0, nstep*nwalker*nvalue*sizeof(double));
  memset(recs_.lnp[1], 0, nstep*nwalker*sizeof(double));
  memset(recs_.newstate[1], 0, nstep*nwalker*sizeof(int));
}

void gather_probability(SentinelQ<Inversion*> &fitq)
{
  auto q = fitq.getNext();
  // find the last node
  while (q->getNext() != nullptr) {
    q = q->getNext();
  }
  auto qlast = q;

  // replace the log probability
  q = fitq.getNext();
  while (q != qlast) {
    Inversion *pfit = q->getData();

    int cur = pfit->recs_.cur;
    int nwalker = pfit->recs_.nwalker;

    for (int k = 0; k < nwalker; ++k) {
      pfit->recs_.lnp[cur][k] = qlast->getData()->recs_.lnp[cur][k];
    }
    q = q->getNext();
  }
}
