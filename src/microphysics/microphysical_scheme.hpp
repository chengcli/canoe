#ifndef SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// athenapp
#include <athena/athena.hpp>

// snap
#include <air_parcel.hpp>

// utils
#include <utils/ndarray.hpp>
#include <utils/parameter_map.hpp>

class MeshBlock;
class ParameterInput;

class MicrophysicalSchemeBase {
 public:
  MicrophysicalSchemeBase(std::string name) : name_(name) {}

  virtual ~MicrophysicalSchemeBase() {}

  virtual void AssembleReactionMatrix(Real *rate, Real **jac,
                                      AirParcel const &air, Real time) = 0;

  virtual void EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt) = 0;

  Real *GetRatePtr() { return rate_; }

  Real **GetJacobianPtr() { return jacobian_; }

  std::string GetName() { return name_; }

 private:
  std::string name_;
  Real *rate_;
  Real **jacobian_;
}

using MicrophysicalSchemePtr = std::shared_ptr<MicrophysicalSchemeBase>;

template <int D>
class MicrophysicalScheme : public MicrophysicalSchemeBase {
 public:
  enum { Size = D };

  MicrophysicalScheme(std::string name, YAML::Node &node)
      : MicrophysicalSchemeBase(name) {
    NewCArray(rate_, Size);
    NewCArray(jacobian_, Size, Size);
  }

  ~MicrophysicalScheme() {
    FreeCArray(rate_);
    FreeCArray(jacobian_);
  }

  protectecd :
      //! parameters
      ParameterMap params_;

  //! indices of species
  std::array<int, D> species_index_;
};

class Kessler94 : public MicrophysicalScheme<3> {
 public:
  using Base = MicrophysicalScheme<3>;

  Kessler94(std::string name, YAML::Node &node);

  ~Kessler94();

  void AssembleReactionMatrix(Real *rate, Real **jac, AirParcel const &air,
                              Real time) override;

  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

 protected:
  ChemistrySolve<Base::Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
