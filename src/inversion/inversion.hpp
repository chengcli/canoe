#ifndef INVERSION_HPP
#define INVERSION_HPP

// C/C++ headers
#include <string>

// Eigen headers
#include <Eigen/Core>

// harp2 headers
#include <configure.hpp>
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class BlockIndex;
class Hydro;
class Thermodynamics;
class Radiation;
template<typename T>
class SentinelQ;

class Inversion {
  friend void gather_probability(SentinelQ<Inversion*> &fitq);
public:
  // functions
  Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name);

  virtual ~Inversion();

	virtual Real LogPosteriorProbability(Radiation *prad, Hydro *phydro,
    Real const *par, Real *val, int k) const
  {
    return 0.;
  }

  virtual void CalculateFitTarget(Radiation const *prad,
    Real *val, int nvalue, int k, int j) const
  {}

  virtual void InitializePositions()
  {}

  virtual void UpdateHydro(Hydro *phydro, ParameterInput *pin) const
  {}

  virtual int getX2Span() const
  {
    return 0.;
  }

  // MCMC functions
	void InitializeChain(int nstep, int nwalker, int ndim, int nvalue);

  void MakeMCMCOutputs(std::string fname);

  void MCMCInit(Radiation *prad, Hydro *phydro);

  void MCMCMove(Radiation *prad, Hydro *phydro);

  void MCMCSave(Hydro *phydro);

  void ResetChain();

  // access functions
  int getDims() const
  {
    return recs_.ndim;
  }

  int getValues() const
  {
    return recs_.nvalue;
  }

  int getWalkers() const
  {
    return recs_.nwalker;
  }

  int getSteps() const
  {
    return recs_.nstep;
  }

  std::string getName() const
  {
    return name_;
  }

  void setX2Indices(int j)
  {
    jl_ = j;
    ju_ = jl_ + getX2Span() - 1;
  }

  Inversion* use(BlockIndex const *p)
  {
    pblock_ = p;
    return this;
  }

  Inversion* use(Thermodynamics const *p)
  {
    pthermo_ = p;
    return this;
  }

protected:
  // name of the inversion
  std::string             name_;

  // fit data
  Eigen::VectorXd         target_;
  Eigen::MatrixXd         icov_;

  // connection
  Coordinates     const*  pcoord_;
  BlockIndex      const*  pblock_;
  Thermodynamics  const*  pthermo_;

	// mcmc initial positions
  Real**                  init_pos_;

  // whether to fit differential observation
  bool                    fit_differential_;

  // j-index range
  int                     jl_, ju_;

private:
  // mcmc variables
  mcmc_opts               opts_;
  mcmc_recs               recs_;
  bool                    mcmc_initialized_;

  // scratch arrays
  Real                    *zz_, *par_;
};

void read_observation_file(
  Eigen::VectorXd &target,
  Eigen::MatrixXd &icov,
  std::string fname
  );

void new_inversion_queue(SentinelQ<Inversion*> &fitq,
  MeshBlock *pmb, ParameterInput *pin,
  BlockIndex *pblock, Thermodynamics *pthermo);

void gather_probability(SentinelQ<Inversion*> &fitq);

#endif
