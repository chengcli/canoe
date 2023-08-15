#ifndef SRC_MICROPHYSICS_MICROPHYSICS_HPP_
#define SRC_MICROPHYSICS_MICROPHYSICS_HPP_

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

class MeshBlock;
class ParameterInput;

class MicrophysicsSystemBase {
 public:
  MicrophysicsSystemBase() {}

  virtual ~MicrophysicsSystemBase() {}

  virtual void AssembleReactionMatrix(Real *rate, Real **jac,
                                      AirParcel const &air, Real time) {}

  virtual void EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt) {}

  virtual void SetSedimentationVelocity(int k, int j, int il, int iu) {}

  Real *GetRatePtr() { return rate_; }

  Real **GetJacobianPtr() { return jacobian_; }

 private:
  Real *rate_;
  Real **jacobian_;
}

using MicrophysicsSystemPtr = std::shared_ptr<MicrophysicsSystemBase>;

template <int D>
class MicrophysicsSystem : public MicrophysicsSystemBase {
 public:
  enum { Dimension = D };

  using Matrix = Eigen::Matrix<Real, D, D>;
  using Vecor = Eigen::Matrix<Real, D, 1>;

  MicrophysicsSystem(MeshBlock *pmb, ParameterInput *pin);

  ~MicrophysicsSystem() {}

  virtual void AssembleReactionMatrix(Real *rate, Real **jac,
                                      AirParcel const &air, Real time) override;

  virtual void EvolveOneStep(AirParcel *air, Real time, Real dt) override;

  virtual void SetSedimentationVelocity(int k, int j, int il, int iu) override;

  protectecd : Vector rate_;
  Matrix jac_;

  //! reaction coefficients
  std::map<std::string, Real> coeffs_;

  //! indices of cloud variables
  std::vector<int> cloud_index_;
};

class Microphysics {
 public:
  enum { NCLOUD_HYDRO = 4 };

  AthenaArray<Real> w, u;

  Microphysics(MeshBlock *pmb, ParameterInput *pin);

  ~Microphysics();

  void EvolveSystems(std::vector<AirParcel> &air_column, Real time, Real dt) {
    for (auto &system : systems_)
      for (auto &air : air_column) {
        system->AssembleReactionMatrix(system->GetRatePtr(),
                                       system->GetJacobianPtr(), air, time);
        system->EvolveOneStep(&air, time, dt);
      }
  }

  void AddFrictionalHeating(std::vector<AirParcel> &air_column);

  void AddSedimentationFlux(AthenaArray<Real> &u, int k, int j, int il, int iu);

 protected:
  AthenaArray<Real> vsed_, vsedf_;

  AthenaArray<Real> hydro_;

  std::vector<MicrophysicsSystemPtr> systems_;

  MeshBlock *pmy_block_;
};

#endif  // SRC_MICROPHYSICS_MICROPHYSICS_HPP_
