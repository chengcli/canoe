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
#include <utils/ndarrays.hpp>
#include <utils/parameter_map.hpp>

// microphysics
#include "chemistry_solver.hpp"

class MeshBlock;
class ParameterInput;

class MicrophysicalSchemeBase {
 public:
  MicrophysicalSchemeBase(std::string name)
      : name_(name), rate_(nullptr), jacobian_(nullptr) {}

  virtual ~MicrophysicalSchemeBase() {}

  virtual void AssembleReactionMatrix(Real *rate, Real **jac,
                                      AirParcel const &air, Real time) = 0;

  virtual void SetSedimentationVelocity(AthenaArray<Real> &vsed, int k, int j,
                                        int il, int iu) = 0;

  virtual void EvolveOneStep(AirParcel *Air, Real time, Real dt) = 0;

  Real *GetRatePtr() { return rate_; }

  Real **GetJacobianPtr() { return jacobian_; }

  std::string GetName() { return name_; }

 protected:
  std::string name_;
  Real *rate_;
  Real **jacobian_;
};

using MicrophysicalSchemePtr = std::shared_ptr<MicrophysicalSchemeBase>;

template <int D>
class MicrophysicalScheme : public MicrophysicalSchemeBase {
 public:
  enum { Size = D };
  using MapVector = Eigen::Map<Eigen::Matrix<Real, Size, 1>>;
  using MapMatrix = Eigen::Map<Eigen::Matrix<Real, Size, Size>>;
  using RealVector = Eigen::Matrix<Real, Size, 1>;
  using RealMatrix = Eigen::Matrix<Real, Size, Size>;

  MicrophysicalScheme(std::string name, YAML::Node const &node)
      : MicrophysicalSchemeBase(name) {
    rate_ = new Real[Size];
    NewCArray(jacobian_, Size, Size);
  }

  ~MicrophysicalScheme() {
    delete[] rate_;
    FreeCArray(jacobian_);
  }

 protected:
  //! parameters
  ParameterMap params_;

  //! indices of species
  std::array<int, D> species_index_;
};

class Kessler94 : public MicrophysicalScheme<3> {
 public:
  using Base = MicrophysicalScheme<3>;

  Kessler94(std::string name, YAML::Node const &node);

  ~Kessler94();

  void AssembleReactionMatrix(Real *rate, Real **jac, AirParcel const &air,
                              Real time) override;

  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

  void SetSedimentationVelocity(AthenaArray<Real> &vsed, int k, int j, int il,
                                int iu) override;

 protected:
  ChemistrySolver<Base::Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
