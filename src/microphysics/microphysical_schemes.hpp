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
#include <virtual_groups.hpp>

// microphysics
#include "chemistry_solver.hpp"

class MeshBlock;
class ParameterInput;

//! \brief virtual base class for all microphysical schemes
class MicrophysicalSchemeBase : public NamedGroup {
 public:  /// constructor and destructor
  MicrophysicalSchemeBase(std::string name) : NamedGroup(name) {}
  virtual ~MicrophysicalSchemeBase() {}

 public:  /// functions
  //! \brief Assemble the reaction matrix
  //!
  //! \param [in] air air parcel provides data to populate the reaction matrix
  //! \param [in] time current simulation time
  virtual void AssembleReactionMatrix(AirParcel const &air, Real time) = 0;

  //! \brief Evolve the air parcel one time step
  //!
  //! \param [in,out] air air parcel to be evolved
  //! \param [in] time current simulation time
  //! \param [in] dt time step
  virtual void EvolveOneStep(AirParcel *air, Real time, Real dt) = 0;

 public:  /// inbound functions
  //! \brief Set the sedimentation velocity from the conserved variables
  //!
  //! \param [in,out] vsed sedimentation velocity
  //! \param [in] phydro root-level hydrodynamic class
  virtual void SetVsedFromConserved(AthenaArray<Real> vsed[3],
                                    Hydro const *phydro, int kl, int ku, int jl,
                                    int ju, int il, int iu) = 0;
};

using MicrophysicalSchemePtr = std::shared_ptr<MicrophysicalSchemeBase>;
using AllMicrophysicalSchemes = std::vector<MicrophysicalSchemePtr>;

//! \brief factory class for contructing microphysical schemes
class MicrophysicalSchemesFactory {
 public:
  static AllMicrophysicalSchemes Create(MeshBlock *pmb, ParameterInput *pin);
};

//! \brief base class for all microphysical schemes
template <int D>
class MicrophysicalScheme : public MicrophysicalSchemeBase,
                            public ParameterGroup,
                            public SpeciesIndexGroup {
 public:
  // needed for Eigen small matrix
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  enum { Size = D };
  using RealVector = Eigen::Matrix<Real, D, 1>;
  using RealMatrix = Eigen::Matrix<Real, D, D>;

 public:  /// constructor and destructor
  MicrophysicalScheme(std::string name, YAML::Node const &node)
      : MicrophysicalSchemeBase(name) {
    if (node["parameters"]) SetRealsFrom(node["parameters"]);

    if (node["dependent-species"]) {
      auto species = node["dependent-species"].as<std::vector<std::string>>();
      SetSpeciesIndex(species);
    }
  }

  virtual ~MicrophysicalScheme() {}

 public:  /// member functions
  Real const *GetRatePtr() const { return rate_.data(); }
  Real const *GetJacobianPtr() const { return jacb_.data(); }

 public:  /// inbound functions
  void SetVsedFromConserved(AthenaArray<Real> vsed[3], Hydro const *phydro,
                            int kl, int ku, int jl, int ju, int il,
                            int iu) override {
    for (auto n : GetMyCloudIndices())
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          for (int i = il; i <= iu; ++i) {
            vsed[0](n, k, j, i) = 0.0;
            vsed[1](n, k, j, i) = 0.0;
            vsed[2](n, k, j, i) = 0.0;
          }
  }

 protected:
  //! rate and jacobian
  Eigen::Matrix<Real, D, 1> rate_;
  Eigen::Matrix<Real, D, D> jacb_;
};

//! \brief Kessler (1994) microphysical scheme
class Kessler94 : public MicrophysicalScheme<3> {
 public:
  using Base = MicrophysicalScheme<3>;
  enum { Size = Base::Size };

 public:  /// constructor and destructor
  Kessler94(std::string name, YAML::Node const &node);
  ~Kessler94();

 public:  /// functions
  void AssembleReactionMatrix(AirParcel const &air, Real time) override;
  void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

  /// inbound functions
  void SetVsedFromConserved(AthenaArray<Real> vsed[3], Hydro const *phydro,
                            int kl, int ku, int jl, int ju, int il,
                            int iu) override;

 protected:
  //! chemistry solver
  ChemistrySolver<Size> solver_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICAL_SCHEMES_HPP_
