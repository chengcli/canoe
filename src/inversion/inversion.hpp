#ifndef SRC_INVERSION_INVERSION_HPP_
#define SRC_INVERSION_INVERSION_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// canone
#include <configure.hpp>

// Eigen
#include <Eigen/Core>

// athena
#include <athena/athena.hpp>

// inversion
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class BlockIndex;
class Hydro;
class Thermodynamics;
class Radiation;

class Inversion {
 public:
  static std::vector<std::unique_ptr<Inversion>> NewInversionQueue(
      MeshBlock *pmb, ParameterInput *pin);

  virtual ~Inversion();

  virtual Real LogPosteriorProbability(Radiation *prad, Hydro *phydro,
                                       Real const *par, Real *val,
                                       int k) const {
    return 0.;
  }

  virtual void CalculateFitTarget(Radiation const *prad, Real *val, int nvalue,
                                  int k, int j) const {}

  virtual void InitializePositions() {}

  virtual void UpdateHydro(Hydro *phydro, ParameterInput *pin) const {}

  virtual int getX2Span() const { return 0.; }

  // MCMC functions
  void InitializeChain(int nstep, int nwalker, int ndim, int nvalue);

  void MakeMCMCOutputs(std::string fname);

  void MCMCInit(Radiation *prad, Hydro *phydro);

  void MCMCMove(Radiation *prad, Hydro *phydro);

  void MCMCSave(Hydro *phydro);

  void ResetChain();

  // access functions
  int GetDims() const { return recs_.ndim; }

  int GetValues() const { return recs_.nvalue; }

  int GetWalkers() const { return recs_.nwalker; }

  int GetSteps() const { return recs_.nstep; }

  void SetLogProbability(int k, Real lnp) { recs_.lnp[recs_.cur][k] = lnp; }

  Real GetLogProbability(int k) const { return recs_.lnp[recs_.cur][k]; }

  std::string GetName() const { return name_; }

  void setX2Indices(int j) {
    jl_ = j;
    ju_ = jl_ + getX2Span() - 1;
  }

 protected:
  Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name);

  // name of the inversion
  std::string name_;

  // fit data
  Eigen::VectorXd target_;
  Eigen::MatrixXd icov_;

  MeshBlock *pmy_block_;

  // mcmc initial positions
  Real **init_pos_;

  // whether to fit differential observation
  bool fit_differential_;

  // j-index range
  int jl_, ju_;

 private:
  // mcmc variables
  mcmc_opts opts_;
  mcmc_recs recs_;
  bool mcmc_initialized_;

  // scratch arrays
  Real *zz_, *par_;
};

using InversionPtr = std::unique_ptr<Inversion>;
using InversionQueue = std::vector<std::unique_ptr<Inversion>>;

#endif  //  SRC_INVERSION_INVERSION_HPP_
