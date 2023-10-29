#ifndef SRC_SNAP_TURBULENCE_TURBULENCE_MODEL_HPP_
#define SRC_SNAP_TURBULENCE_TURBULENCE_MODEL_HPP_

// C/C++
#include <memory>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/bvals/cc/bvals_cc.hpp>

class MeshBlock;
class ParameterInput;

//! \class TurbulenceModel
//  \brief

class TurbulenceModel {
 public:
  TurbulenceModel(MeshBlock *pmb, ParameterInput *pin);
  virtual ~TurbulenceModel();

  // access members
  AthenaArray<Real> w, u;

  // public data
  AthenaArray<Real> mut;  // dynamic turbulent viscosity

  // public functions:
  virtual void DriveTurbulence(Real dt) {}
  virtual void Initialize() {}
  virtual void SetDiffusivity(AthenaArray<Real> &nu, AthenaArray<Real> &kappa,
                              const AthenaArray<Real> &w,
                              const AthenaArray<Real> &bc, int il, int iu,
                              int jl, int ju, int kl, int ku) {}

 protected:
  MeshBlock *pmy_block;
};

class KEpsilonTurbulence : public TurbulenceModel {
 public:
  KEpsilonTurbulence(MeshBlock *pmb, ParameterInput *pin);
  ~KEpsilonTurbulence() {}

  void DriveTurbulence(Real dt);
  void Initialize();
  void SetDiffusivity(AthenaArray<Real> &nu, AthenaArray<Real> &kappa,
                      const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
                      int il, int iu, int jl, int ju, int kl, int ku);

 private:
  Real cmu_, c1_, c2_, sigk_, sige_;
};

using TurbulenceModelPtr = std::shared_ptr<TurbulenceModel>;

class TurbulenceFactory {
 public:
  static TurbulenceModelPtr Create(MeshBlock *pmb, ParameterInput *pin);
};

#endif  // SRC_SNAP_TURBULENCE_TURBULENCE_MODEL_HPP_
