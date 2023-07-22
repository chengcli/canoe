#ifndef TURBULENCE_MODEL_HPP_
#define TURBULENCE_MODEL_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file scalars.hpp
//  \brief definitions for TurbulenceModel class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;

//! \class TurbulenceModel
//  \brief

class TurbulenceModel {
 public:
  TurbulenceModel(MeshBlock *pmb, ParameterInput *pin, int nvar = 1);
  virtual ~TurbulenceModel() {}

  // public data:
  // "conserved vars" = extensive variable
  AthenaArray<Real> s, s1, s2;  // (no more than MAX_NREGISTER allowed)
  // "primitive vars" = intensive variable
  AthenaArray<Real> r;          // , r1;
  AthenaArray<Real> s_flux[3];  // face-averaged flux vector
  AthenaArray<Real> mut;        // dynamic turbulent viscosity

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_s_, coarse_r_;
  int refinement_idx{-1};

  CellCenteredBoundaryVariable sbvar;

  // public functions:
  void addFluxDivergence(const Real wght, AthenaArray<Real> &s_out);
  void calculateFluxes(AthenaArray<Real> &s, const int order);
  void ConservedToPrimitive(AthenaArray<Real> &s, AthenaArray<Real> const &w,
                            AthenaArray<Real> &r, Coordinates *pco, int il,
                            int iu, int jl, int ju, int kl, int ku);
  void PrimitiveToConserved(AthenaArray<Real> &r, AthenaArray<Real> const &w,
                            AthenaArray<Real> &s, Coordinates *pco, int il,
                            int iu, int jl, int ju, int kl, int ku);
  void computeUpwindFlux(const int k, const int j, const int il,
                         const int iu,  // CoordinateDirection dir,
                         AthenaArray<Real> &rl, AthenaArray<Real> &rr,
                         AthenaArray<Real> &mass_flx,
                         AthenaArray<Real> &flx_out);
  void applyBoundaryCondition(AthenaArray<Real> &r, AthenaArray<Real> &s,
                              AthenaArray<Real> const &w, Coordinates *pco);

  virtual void driveTurbulence(AthenaArray<Real> &s, AthenaArray<Real> const &r,
                               AthenaArray<Real> const &w, Real dt) {}
  virtual void Initialize(AthenaArray<Real> const &w) {}

  virtual void setDiffusivity(AthenaArray<Real> &nu, AthenaArray<Real> &kappa,
                              const AthenaArray<Real> &w,
                              const AthenaArray<Real> &bc, int il, int iu,
                              int jl, int ju, int kl, int ku) {}

  // restart functions
  size_t getRestartDataSizeInBytes();
  size_t dumpRestartData(char *pdst);
  size_t loadRestartData(char *psrc);

 protected:
  MeshBlock *pmy_block;
  // scratch space used to compute fluxes
  // 2D scratch arrays
  AthenaArray<Real> rl_, rr_, rlb_;
  // 1D scratch arrays
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;
};

class KEpsilonTurbulence : public TurbulenceModel {
 public:
  KEpsilonTurbulence(MeshBlock *pmb, ParameterInput *pin);
  ~KEpsilonTurbulence() {}
  void driveTurbulence(AthenaArray<Real> &s, AthenaArray<Real> const &r,
                       AthenaArray<Real> const &w, Real dt);
  void Initialize(AthenaArray<Real> const &w);

  void setDiffusivity(AthenaArray<Real> &nu, AthenaArray<Real> &kappa,
                      const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
                      int il, int iu, int jl, int ju, int kl, int ku);

 private:
  Real cmu_, c1_, c2_, sigk_, sige_;
};

#endif  // TURBULENCE_MODEL_HPP_
