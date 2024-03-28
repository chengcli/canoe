#ifndef SRC_INVERSION_INVERSION_HPP_
#define SRC_INVERSION_INVERSION_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <configure.hpp>
#include <virtual_groups.hpp>

// inversion
#include "mcmc.hpp"

class MeshBlock;
class ParameterInput;

//! \brief Base class for inversion
//!
//! This class is the base class for inversion.
//! It provides the interface for one function, UpdateModel.
//! The UpdateModel function takes a vector of parameters and updates the model.
//! The parameters are to be adjusted for a better fit to the data.
//! The UpdateModel function may take several steps to update the model.
//! The results of each step are stored from index jl_ to ju_ (exclusive).
class Inversion : public ParameterGroup {
 public:
  enum {
    Steps = 0,
  }

  /// Constructor and destructor
  Inversion(MeshBlock *pmb, ParameterInput *pin);
  virtual ~Inversion();

  virtual void UpdateModel(std::vector<Real> const &par) const {}

  void SetStepRange(int js, int je) {
    jl_ = js;
    ju_ = je;
  }

 protected:
  // step index range (j-direction)
  int jl_, ju_;

  //! pointer to parent MeshBlock
  MeshBlock const *pmy_block_;
};

class CompositionInversion : public Inversion {
 public:
  CompositionInversion(MeshBlock *pmb, ParameterInput *pin);
  ~CompositionInversion();

  int GetSteps() const override { return idx_.size(); }

  void UpdateConcentration(Hydro *phydro, Real *Xp, int k, int jl,
                           int ju) const;

 protected:
  // inversion variable ids
  std::vector<int> idx_;

  // prior standard deviation
  Real Xstd_[1 + NVAPOR];
};

class ProfileInversion : public Inversion {
 public:
  ProfileInversion(MeshBlock *pmb, ParameterInput *pin);
  ~ProfileInversion() {}

  size_t samples() const { return plevel_.size() - 2; }

  void UpdateModel(std::vector<Real> const &par) const override;

  void UpdateProfiles(Hydro *phydro, Real **XpSample, int k, int jl,
                      int ju) const;

  void ConvectiveAdjustment(Hydro *phydro, int k, int ju) const;

  int GetSteps() const override { return idx_.size(); }

 protected:
  // pressure levels
  std::vector<Real> plevel_;

  // inversion variable ids
  std::vector<int> idx_;

  // hyper-parameters
  Real chi_;
  Real Xstd_[1 + NVAPOR];
  Real Xlen_[1 + NVAPOR];
};

using InversionPtr = std::shared_ptr<Inversion>;

class InversionsFactory {
 public:
  static std::vector<InversionPtr> CreateFrom(MeshBlock *pmb,
                                              ParameterInput *pin);
};

#endif  //  SRC_INVERSION_INVERSION_HPP_
