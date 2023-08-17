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

class MeshBlock;
class ParameterInput;

class MicrophysicalSchemeBase {
 public:
  MicrophysicalSchemeBase() {}

  virtual ~MicrophysicalSchemeBase() {}

  virtual void AssembleReactionMatrix(Real *rate, Real **jac,
                                      AirParcel const &air, Real time) = 0;

  virtual void EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt) = 0;

  Real *GetRatePtr() { return rate_; }

  Real **GetJacobianPtr() { return jacobian_; }

 private:
  Real *rate_;
  Real **jacobian_;
}

using MicrophysicalSchemePtr = std::shared_ptr<MicrophysicalSchemeBase>;

template <int D>
class MicrophysicalScheme : public MicrophysicalSchemeBase {
 public:
  enum { Size = D };

  MicrophysicalScheme(YAML::Node &node) {
    NewCArray(rate_, Size);
    NewCArray(jacobian_, Size, Size);
  }

  ~MicrophysicalScheme() {
    FreeCArray(rate_);
    FreeCArray(jacobian_);
  }

  protectecd :
      //! reaction coefficients
      std::map<std::string, Real>
          coeffs_;

  //! indices of species
  std::array<int, D> species_index_;
};

class Kessler94 : public MicrophysicalScheme<3> {
 public:
  Kessler94(YAML::Node &node);

  ~Kessler94();

  void AssembleReactionMatrix(Real *rate, Real **jac, AirParcel const &air,
                              Real time) override;

  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

 protected:
  ChemistrySolve<3> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
