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

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
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

  virtual void SetSedimentationVelocityX1(AthenaArray<Real> &vsed, int kl,
                                          int ku, int jl, int ju, int il,
                                          int iu) = 0;

  virtual void SetSedimentationVelocityX2(AthenaArray<Real> &vsed, int kl,
                                          int ku, int jl, int ju, int il,
                                          int iu) = 0;

  virtual void SetSedimentationVelocityX3(AthenaArray<Real> &vsed, int kl,
                                          int ku, int jl, int ju, int il,
                                          int iu) = 0;

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

  void SetSedimentationVelocityX1(AthenaArray<Real> &vsed,
                                  Meshblock *pmb) override {
    int ks = pmb->ks, js = pmb->js, is = pmb->is;
    int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

    for (auto n : species_index_)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) vsed(n - NHYDRO, k, j, i) = 0.0;
  }

  void SetSedimentationVelocityX2(AthenaArray<Real> &vsed,
                                  MeshBlock *pmb) override {
    int ks = pmb->ks, js = pmb->js, is = pmb->is;
    int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

    for (auto n : species_index_)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) vsed(n - NHYDRO, k, j, i) = 0.0;
  }

  void SetSedimentationVelocityX3(AthenaArray<Real> &vsed,
                                  MeshBlock *pmb) override {
    int ks = pmb->ks, js = pmb->js, is = pmb->is;
    int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

    for (auto n : species_index_)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) vsed(n - NHYDRO, k, j, i) = 0.0;
  }

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

  void SetSedimentationVelocityX1(AthenaArray<Real> &vsed,
                                  MeshBlock *pmb) override;

 protected:
  ChemistrySolver<Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEME_HPP_
