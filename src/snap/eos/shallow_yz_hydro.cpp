//! \file shallow_yz.cpp
//  \brief implements functions in class EquationOfState for shallow water
//  hydrodynamics

// C/C++
#include <cmath>  // sqrt()

// athena
#include <athena/athena.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <athena/stride_iterator.hpp>

// canoe
#include <configure.hpp>

// exo3
#include <exo3/cubed_sphere_utility.hpp>
#include <exo3/gnomonic_equiangle.hpp>

namespace cs = CubedSphereUtility;

// EquationOfState constructor
//
EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block_(pmb),
      density_floor_{
          pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024 * float_min))},
      pressure_floor_{
          pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024 * float_min))} {}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku)
// \brief Converts conserved into primitive variables in shallow water.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old,
    const FaceField &b, AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real &ud = cons(IDN, k, j, i);

        // apply density floor, without changing momentum or energy
        ud = (ud > density_floor_) ? ud : density_floor_;
        prim(IDN, k, j, i) = ud;
        prim(IVX, k, j, i) = 0.;
        prim(IVY, k, j, i) = cons(IVY, k, j, i);
        prim(IVZ, k, j, i) = cons(IVZ, k, j, i);
        Real di = 1. / ud;

#ifdef CUBED_SPHERE
        cs::CovariantToContravariant(
            prim.at(k, j, i),
            static_cast<GnomonicEquiangle *>(pco)->GetCosineCell(k, j));
#endif
        prim(IVY, k, j, i) *= di;
        prim(IVZ, k, j, i) *= di;
      }
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real>
// &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates
//           *pco, int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
                                           const AthenaArray<Real> &bc,
                                           AthenaArray<Real> &cons,
                                           Coordinates *pco, int il, int iu,
                                           int jl, int ju, int kl, int ku) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        const Real &wd = prim(IDN, k, j, i);
        cons(IDN, k, j, i) = wd;
        cons(IVX, k, j, i) = 0.;
        cons(IVY, k, j, i) = prim(IVY, k, j, i);
        cons(IVZ, k, j, i) = prim(IVZ, k, j, i);

#ifdef CUBED_SPHERE
        cs::ContravariantToCovariant(
            cons.at(k, j, i),
            static_cast<GnomonicEquiangle *>(pco)->GetCosineCell(k, j));
#endif
        cons(IVY, k, j, i) *= wd;
        cons(IVZ, k, j, i) *= wd;
      }
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns shallow water gravity wave speed given vector of primitive
// variables

Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  return std::sqrt(prim[IDN]);
}

//---------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim,
//           int k, int j, int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface
// states

void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k,
                                           int j, int i) {
  Real &w_d = prim(IDN, i);

  // apply density floor
  w_d = (w_d > density_floor_) ? w_d : density_floor_;

  return;
}

void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real> &prim,
                                                    AthenaArray<Real> &cons,
                                                    AthenaArray<Real> &bcc,
                                                    int k, int j, int i) {
  Real &w_d = prim(IDN, k, j, i);
  Real &u_d = cons(IDN, k, j, i);

  // apply density floor
  w_d = (w_d > density_floor_) ? w_d : density_floor_;
  u_d = w_d;
}
