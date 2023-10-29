#ifndef SRC_INVERSION_INVERSION_HPP_
#define SRC_INVERSION_INVERSION_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// Eigen
#include <Eigen/Core>

// athena
#include <athena/athena.hpp>

// canoe
#include <configure.hpp>
#include <virtual_groups.hpp>

// inversion
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Hydro;
class Thermodynamics;
class Radiation;

class Inversion : public NamedGroup,
                  //public RestartGroup,
                  public FITSOutputGroup {
 public:
  /// Constructor and destructor
  Inversion(MeshBlock *pmb, ParameterInput *pin, std::string name);
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

  void setX2Indices(int j) {
    jl_ = j;
    ju_ = jl_ + getX2Span() - 1;
  }

  /// MeshOutputGroup functions
  bool ShouldFITSOutput(std::string variable_name) const override { return true;}
  void LoadFITSOutputData(OutputType *pod, int *num_vars) const override {}

 protected:
  // name of the inversion
  std::string name_;

  // fit data
  Eigen::VectorXd target_;
  Eigen::MatrixXd icov_;

  // mcmc initial positions
  Real **init_pos_;

  // whether to fit differential observation
  bool fit_differential_;

  // j-index range
  int jl_, ju_;

  //! pointer to parent MeshBlock
  MeshBlock const *pmy_block_;

 private:
  // mcmc variables
  mcmc_opts opts_;
  mcmc_recs recs_;
  bool mcmc_initialized_;

  // scratch arrays
  Real *zz_, *par_;
};

using InversionPtr = std::shared_ptr<Inversion>;
using AllInversions = std::vector<std::shared_ptr<Inversion>>;

class InversionsFactory {
 public:
  static AllInversions CreateAllInversions(MeshBlock *pmb, ParameterInput *pin);
};

namespace InversionHelper {

void read_observation_file(Eigen::VectorXd *target, Eigen::MatrixXd *icov,
                           std::string fname);
void gather_probability(std::vector<Inversion *> const &fitq);

}  // namespace InversionHelper

#endif  //  SRC_INVERSION_INVERSION_HPP_
