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
#include <configure.hpp>

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
        Real di = 1. / ud;

        const Real u = 0.;
        const Real v = cons(IVY, k, j, i) * di;
        const Real w = cons(IVZ, k, j, i) * di;

#if defined(CUBED_SPHERE)
        Real x = tan(pco->x2v(j));
        Real y = tan(pco->x3v(k));
        Real C = sqrt(1. + x * x);
        Real D = sqrt(1. + y * y);
        Real cth = -x * y / (C * D);
        Real sth2 = 1. - cth * cth;

        // Raise to contravarient velocity
        prim(IVY, k, j, i) = v / sth2 - w * cth / sth2;
        prim(IVZ, k, j, i) = -v * cth / sth2 + w / sth2;
#else
        prim(IVY, k, j, i) = v;
        prim(IVZ, k, j, i) = w;
#endif
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

        const Real u = 0.;
        const Real v = prim(IVY, k, j, i) * wd;
        const Real w = prim(IVZ, k, j, i) * wd;

#if defined(CUBED_SPHERE)
        Real x = tan(pco->x2v(j));
        Real y = tan(pco->x3v(k));
        Real C = sqrt(1. + x * x);
        Real D = sqrt(1. + y * y);
        Real cth = -x * y / (C * D);

        cons(IVY, k, j, i) = v + w * cth;
        cons(IVZ, k, j, i) = w + v * cth;

#else
        cons(IVY, k, j, i) = v;
        cons(IVZ, k, j, i) = w;
#endif
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
