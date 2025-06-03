// C/C++
#include <cmath>  // sqrt()
#include <iomanip>
#include <sstream>
#include <stdexcept>

// application
#include <application/exceptions.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.h>

#include <impl.hpp>

// exo3
#include <exo3/exo3.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

#include "eos_helper.hpp"

// checks
#include <checks.hpp>

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock* pmb, ParameterInput* pin)
    : pmy_block_(pmb),
      gamma_{pin->GetReal("hydro", "gamma")},
      density_floor_{
          pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024 * float_min))},
      pressure_floor_{
          pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024 * float_min))},
      scalar_floor_{pin->GetOrAddReal("hydro", "sfloor", 0.)} {}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real>& cons, const AthenaArray<Real>& prim_old,
    const FaceField& b, AthenaArray<Real>& prim, AthenaArray<Real>& bcc,
    Coordinates* pco, int il, int iu, int jl, int ju, int kl, int ku) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;

  // apply_vapor_limiter(&cons, pmy_block_);

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real& u_d = cons(IDN, k, j, i);
        Real& u_m1 = cons(IM1, k, j, i);
        Real& u_m2 = cons(IM2, k, j, i);
        Real& u_m3 = cons(IM3, k, j, i);
        Real& u_e = cons(IEN, k, j, i);

        Real& w_d = prim(IDN, k, j, i);
        Real& w_vx = prim(IVX, k, j, i);
        Real& w_vy = prim(IVY, k, j, i);
        Real& w_vz = prim(IVZ, k, j, i);
        Real& w_p = prim(IPR, k, j, i);

        Real density = 0.;
        for (int n = 0; n < IVX; ++n) {
          density += cons(n, k, j, i);
        }
        // total density
        w_d = std::max(density_floor_, density);
        Real di = 1. / density;

        // mass mixing ratio
        for (int n = 1; n < IVX; ++n)
          prim(n, k, j, i) = std::max(scalar_floor_, cons(n, k, j, i) * di);

        w_vx = u_m1 * di;
        w_vy = u_m2 * di;
        w_vz = u_m3 * di;

        // covariant to contravariant
        vec_raise_inplace(prim.at(k, j, i), pco->m.at(k, j, i));

        // internal energy
        Real KE = 0.5 * (u_m1 * w_vx + u_m2 * w_vy + u_m3 * w_vz);
        w_p = pthermo->IntEngToPres(prim.at(k, j, i), u_e - KE);

        // apply pressure floor, correct total energy
        u_e =
            (w_p > pressure_floor_)
                ? u_e
                : pthermo->PresToIntEng(prim.at(k, j, i), pressure_floor_) + KE;
        w_p = (w_p > pressure_floor_) ? w_p : pressure_floor_;
      }
    }
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real>
// &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates
//           *pco, int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(const AthenaArray<Real>& prim,
                                           const AthenaArray<Real>& bc,
                                           AthenaArray<Real>& cons,
                                           Coordinates* pco, int il, int iu,
                                           int jl, int ju, int kl, int ku) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real& u_d = cons(IDN, k, j, i);
        Real& u_m1 = cons(IM1, k, j, i);
        Real& u_m2 = cons(IM2, k, j, i);
        Real& u_m3 = cons(IM3, k, j, i);
        Real& u_e = cons(IEN, k, j, i);

        const Real& w_d = prim(IDN, k, j, i);
        const Real& w_vx = prim(IVX, k, j, i);
        const Real& w_vy = prim(IVY, k, j, i);
        const Real& w_vz = prim(IVZ, k, j, i);
        const Real& w_p = prim(IPR, k, j, i);

        // density
        u_d = w_d;
        for (int n = 1; n < IVX; ++n) {
          cons(n, k, j, i) = prim(n, k, j, i) * w_d;
          cons(IDN, k, j, i) -= cons(n, k, j, i);
        }

        // momentum
        u_m1 = w_vx * w_d;
        u_m2 = w_vy * w_d;
        u_m3 = w_vz * w_d;

        // contravariant to covariant
        vec_lower_inplace(cons.at(k, j, i), pco->m.at(k, j, i));

        // total energy
        Real KE = 0.5 * (u_m1 * w_vx + u_m2 * w_vy + u_m3 * w_vz);
        u_e = pthermo->PresToIntEng(prim.at(k, j, i), w_p) + KE;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  auto pthermo = Thermodynamics::GetInstance();
  Real gamma = pthermo->GetGamma(prim);
  return std::sqrt(gamma * prim[IPR] / prim[IDN]);
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int
// k, int j,
//                                                 =int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface
// states
void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real>& prim, int k,
                                           int j, int i) {
  Real& w_d = prim(IDN, i);
  Real& w_p = prim(IPR, i);

  // apply (prim) density floor
  w_d = (w_d > density_floor_) ? w_d : density_floor_;
  // apply composition floors
  for (int n = 1; n < IVX; ++n) {
    prim(n, i) = std::max(scalar_floor_, prim(n, i));
  }
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ? w_p : pressure_floor_;

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real>
// &prim,
//           AthenaArray<Real> &cons, FaceField &b, int k, int j, int i) {
// \brief Apply pressure (prim) floor and correct energy (cons) (typically after
// W(U))
void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real>& prim,
                                                    AthenaArray<Real>& cons,
                                                    AthenaArray<Real>& bcc,
                                                    int k, int j, int i) {
  auto pthermo = Thermodynamics::GetInstance();
  Real& w_d = prim(IDN, k, j, i);
  Real& w_p = prim(IPR, k, j, i);

  Real& u_d = cons(IDN, k, j, i);
  Real& u_e = cons(IEN, k, j, i);

  // apply (prim) density floor, without changing momentum or energy
  w_d = (w_d > density_floor_) ? w_d : density_floor_;
  // ensure cons density matches
  u_d = w_d;

  // apply composition floors
  for (int n = 1; n < IVX; ++n) {
    prim(n, k, j, i) = std::max(scalar_floor_, prim(n, k, j, i));
    cons(n, k, j, i) = prim(n, k, j, i) * w_d;
    u_d -= cons(n, k, j, i);
  }

  auto vl = vec_lower(prim.at(k, j, i), pmy_block_->pcoord->m.at(k, j, i));

  Real e_k = 0.5 * w_d *
             (prim(IVX, k, j, i) * vl[IVX] + prim(IVY, k, j, i) * vl[IVY] +
              prim(IVZ, k, j, i) * vl[IVZ]);

  // apply pressure floor, correct total energy
  u_e = (w_p > pressure_floor_)
            ? u_e
            : pthermo->PresToIntEng(prim.at(k, j, i), pressure_floor_) + e_k;
  w_p = (w_p > pressure_floor_) ? w_p : pressure_floor_;

  return;
}
