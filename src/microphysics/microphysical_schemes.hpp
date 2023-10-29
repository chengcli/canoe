#ifndef SRC_MICROPHYSICS_MICROPHYSICAL_SCHEMES_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICAL_SCHEMES_HPP_

// C/C++
#include <map>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>
#include <Eigen/Dense>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <index_map.hpp>
#include <virtual_groups.hpp>

// utils
#include <utils/parameter_map.hpp>

// microphysics
#include "chemistry_solver.hpp"

class MeshBlock;
class ParameterInput;

class MicrophysicalSchemeBase : public NamedGroup {
 public:
  /// constructor and destructor
  MicrophysicalSchemeBase(std::string name);
  virtual ~MicrophysicalSchemeBase() {}

  /// functions
  virtual void AssembleReactionMatrix(AirParcel const &air, Real time) = 0;
  virtual void EvolveOneStep(AirParcel *Air, Real time, Real dt) = 0;

  /// inbound functions
  virtual void SetSedimentationVelocityFromConserved(Hydro const *phydro) = 0;

 protected:
  AthenaArray<Real> vsed_shallow_[3];
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

  /// constructor and destructor
  MicrophysicalScheme(std::string name, YAML::Node const &node)
      : MicrophysicalSchemeBase(name) {
    if (node["parameters"]) params_ = ToParameterMap(node["parameters"]);

    std::vector<std::string> species;
    if (node["dependent-species"])
      species = node["dependent-species"].as<std::vector<std::string>>();

    auto pindex = IndexMap::GetInstance();

    for (int i = 0; i < Size; ++i) {
      species_index_[i] = pindex->GetSpeciesId(species[i]);
    }
  }

  virtual ~MicrophysicalScheme() {}

  /// functions
  Real const *GetRatePtr() const { return rate_.data(); }
  Real const *GetJacobianPtr() const { return jacb_.data(); }

  /// inbound functions
  void SetSedimentationVelocityFromConserved(Hydro const *phydro, int kl,
                                             int ku, int jl, int ju, int il,
                                             int iu) override {
    for (auto n : species_index_)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          for (int i = il; i <= iu; ++i) {
            vsed_shallow_[0](n - NHYDRO, k, j, i) = 0.0;
            vsed_shallow_[1](n - NHYDRO, k, j, i) = 0.0;
            vsed_shallow_[2](n - NHYDRO, k, j, i) = 0.0;
          }
  }

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

  /// constructor and destructor
  Kessler94(std::string name, YAML::Node const &node);
  ~Kessler94();

  /// functions
  void AssembleReactionMatrix(AirParcel const &air, Real time) override;
  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

  /// inbound functions
  void SetSedimentationVelocityFromConservedX1(Hydro const *phydro, int kl,
                                               int ku, int jl, int ju, int kl,
                                               int ku) override;

 protected:
  ChemistrySolver<Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEMES_HPP_
