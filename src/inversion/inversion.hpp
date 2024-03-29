#ifndef SRC_INVERSION_INVERSION_HPP_
#define SRC_INVERSION_INVERSION_HPP_

// C/C++
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <virtual_groups.hpp>

class MeshBlock;

//! \brief Base class for inversion
//!
//! This class is the base class for inversion.
//! It provides the interface for two functions, UpdateModel and GetSteps.
//! The UpdateModel function takes a vector of parameters and updates the model.
//! The parameters are to be adjusted for a better fit to the data.
//! The UpdateModel function may take several steps to update the model.
//! The results of each step are stored from index jl_ to ju_ (exclusive).
class Inversion : public NamedGroup,
                  public ParameterGroup,
                  public SpeciesIndexGroup {
 public:
  /// Constructor and destructor
  Inversion(std::string name) : NamedGroup(name), jl_(0), ju_(0) {}

  virtual ~Inversion() {}

  virtual void UpdateModel(MeshBlock *pmb, std::vector<Real> const &par,
                           int k) const {}

  virtual int GetSteps() const { return 0; }

  void SetStepRange(int js, int je) {
    jl_ = js;
    ju_ = je;
  }

 protected:
  // step index range (j-direction)
  int jl_, ju_;
};

class CompositionInversion : public Inversion {
 public:
  CompositionInversion(YAML::Node const &node);
  ~CompositionInversion() {}

  int GetSteps() const override { return GetMySpeciesIndices().size(); }

  void UpdateModel(MeshBlock *pmb, std::vector<Real> const &par,
                   int k) const override;

 protected:
  // prior standard deviation
  Real Xstd_[1 + NVAPOR];
};

class ProfileInversion : public Inversion {
 public:
  ProfileInversion(YAML::Node const &node);
  ~ProfileInversion() {}

  void UpdateModel(MeshBlock *pmb, std::vector<Real> const &par,
                   int k) const override;

  void UpdateProfiles(Hydro *phydro, Real **XpSample, int k, int jl,
                      int ju) const;

  int GetSteps() const override { return GetMySpeciesIndices().size() + 1; }

 protected:
  void enforceStability(AirColumn &ac) const;

 protected:
  // pressure levels
  std::vector<Real> plevel_;

  // hyper-parameters
  Real chi_;
  std::vector<Real> Xstd_;
  std::vector<Real> Xlen_;
};

using InversionPtr = std::shared_ptr<Inversion>;

class InversionFactory {
 public:
  static std::vector<InversionPtr> CreateFrom(YAML::Node const &node);
};

#endif  //  SRC_INVERSION_INVERSION_HPP_
