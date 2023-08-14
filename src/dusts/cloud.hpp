#ifndef SRC_DUSTS_CLOUD_HPP_
#define SRC_DUSTS_CLOUD_HPP_

// C/C++
#include <map>
#include <memory>
#include <vector>

// athenapp
#include <athena/athena.hpp>

class MeshBlock;
class ParameterInput;

template <int D>
class MicrophysicsSystem {
 public:
  enum { Dimension = D };

  using Matrix = Eigen::Matrix<Real, D, D>;
  using Vecor = Eigen::Matrix<Real, D, 1>;

  MicrophysicsSystem(MeshBlock *pmb, ParameterInput *pin);
  virtual ~MicrophysicsSystem() {}

  void AssembleReactionMatrix(AirParcel const &air, Real time);

  protectecd : Vector rate_;
  Matrix jac_;

  //! reaction coefficients
  std::map<std::string, Real> coeffs_;
  std::vector<int> index_;
};

using MicrophysicsSystemPtr = std::shared_ptr<Microphysics>;

class Microphysics {
 public:
  AthenaArray<Real> w, u;

  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();

  void EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt);

 protected:
  std::vector<MicrophysicsSystemPtr> systems_;
  MeshBlock *pmy_block_;
};

#endif  // SRC_CLOUDS_CLOUD_HPP_
