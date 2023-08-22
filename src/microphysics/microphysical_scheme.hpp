#ifndef SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>
#include <Eigen/Dense>

// athenapp
#include <athena/athena.hpp>

// snap
#include <air_parcel.hpp>

// utils
#include <utils/parameter_map.hpp>

// microphysics
#include "chemistry_solver.hpp"

class MeshBlock;
class ParameterInput;

class MicrophysicalSchemeBase {
 public:
  MicrophysicalSchemeBase(std::string name) : name_(name) {}

  virtual ~MicrophysicalSchemeBase() {}

  virtual void AssembleReactionMatrix(AirParcel const &air, Real time) = 0;

  virtual void SetSedimentationVelocity(AthenaArray<Real> &vsed, int k, int j,
                                        int il, int iu) = 0;

  virtual void EvolveOneStep(AirParcel *Air, Real time, Real dt) = 0;

  std::string GetName() { return name_; }

 protected:
  std::string name_;
};

using MicrophysicalSchemePtr = std::shared_ptr<MicrophysicalSchemeBase>;

template <int D>
class MicrophysicalScheme : public MicrophysicalSchemeBase {
 public:
  // needed for Eigen small matrix
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  enum { Size = D };
  using RealVector = Eigen::Matrix<Real, D, 1>;
  using RealMatrix = Eigen::Matrix<Real, D, D>;

  MicrophysicalScheme(std::string name, YAML::Node const &node)
      : MicrophysicalSchemeBase(name) {}

  ~MicrophysicalScheme() {}

  Real const *GetRatePtr() const { return rate_.data(); }

  Real const *GetJacobianPtr() const { return jacb_.data(); }

 protected:
  //! parameters
  ParameterMap params_;

  //! indices of species
  std::array<int, D> species_index_;

  //! rate and jacobian
  Eigen::Matrix<Real, D, 1> rate_;
  Eigen::Matrix<Real, D, D> jacb_;
};

class Kessler94 : public MicrophysicalScheme<3> {
 public:
  using Base = MicrophysicalScheme<3>;
  enum { Size = Base::Size };

  Kessler94(std::string name, YAML::Node const &node);

  ~Kessler94();

  void AssembleReactionMatrix(AirParcel const &air, Real time) override;

  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

  void SetSedimentationVelocity(AthenaArray<Real> &vsed, int k, int j, int il,
                                int iu) override;

 protected:
  ChemistrySolver<Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
