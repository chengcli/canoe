// C++ headers
#include <cmath>  // sqrt()
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <debugger/debugger.hpp>

#include "../meshblock_impl.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "eos_helper.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock* pmb, ParameterInput* pin)
    : pmy_block_(pmb),
      gamma_{pin->GetReal("hydro", "gamma")},
      density_floor_{
          pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024 * float_min))},
      pressure_floor_{
          pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024 * float_min))},
      scalar_floor_{
          pin->GetOrAddReal("hydro", "sfloor", std::sqrt(1024 * float_min))} {}

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
  Real gm1 = GetGamma() - 1.0;
  std::stringstream msg;
  Thermodynamics* pthermo = pmy_block_->pimpl->pthermo;

  apply_vapor_limiter(&cons, pmy_block_);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
#pragma omp simd
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

        // apply density floor, without changing momentum or energy
        u_d = (u_d > density_floor_) ? u_d : density_floor_;

        Real density = 0.;
        for (int n = 0; n <= NVAPOR; ++n) density += cons(n, k, j, i);
        w_d = density;
        Real di = 1. / density;

        // mass mixing ratio
        for (int n = 1; n <= NVAPOR; ++n)
          prim(n, k, j, i) = cons(n, k, j, i) * di;

#if DEBUG_LEVEL > 3
        if (std::isnan(w_d) || (w_d < density_floor_)) {  // IDN may be NAN
          msg << "### FATAL ERROR in function ConservedToPrimitive" << std::endl
              << "Density reaches lowest value: " << w_d << std::endl
              << "At position (" << k << "," << j << "," << i << ") in rank "
              << Globals::my_rank << std::endl;
          for (int ii = std::max(i - 3, il); ii <= std::min(i + 3, iu); ++ii) {
            msg << "i = " << ii << " ";
            for (int jj = std::max(j - 3, jl); jj <= std::min(j + 3, ju); ++jj)
              msg << cons(IDN, k, jj, ii) << " ";
            msg << std::endl;
          }
          ATHENA_ERROR(msg);
        }
#endif

        // Real di = 1.0/u_d;
        w_vx = u_m1 * di;
        w_vy = u_m2 * di;
        w_vz = u_m3 * di;

        // internal energy
        Real KE = 0.5 * di * (u_m1 * u_m1 + u_m2 * u_m2 + u_m3 * u_m3);
        Real fsig = 1., feps = 1.;
        // vapors
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += prim(n, k, j, i) * (pthermo->GetCvRatio(n) - 1.);
          feps += prim(n, k, j, i) * (1. / pthermo->GetMassRatio(n) - 1.);
        }
        w_p = gm1 * (u_e - KE) * feps / fsig;

        // apply pressure floor, correct total energy
        u_e = (w_p > pressure_floor_)
                  ? u_e
                  : ((pressure_floor_ / gm1) * fsig / feps + KE);
        w_p = (w_p > pressure_floor_) ? w_p : pressure_floor_;

#if DEBUG_LEVEL > 3
        if (std::isnan(w_p) || (w_p < pressure_floor_)) {
          msg << "### FATAL ERROR in function ConservedToPrimitive" << std::endl
              << "Pressure reaches lowest value: " << w_p << std::endl
              << "At position (" << k << "," << j << "," << i << ") in rank "
              << Globals::my_rank << std::endl;
          ATHENA_ERROR(msg);
        }
#endif
      }
    }
  }

  // #if DEBUG_LEVEL > 0
  //   Debugger *pdbg = pmy_block_->pdebug;
  //   pdbg = pdbg->StartTracking("EquationOfStates::ConservedToPrimitive");
  //   pdbg->Track3D("rho", IsPositive, prim, IDN);
  //   pdbg->Track3D("pres", IsPositive, prim, IPR);
  // #endif

  return;
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
  Real igm1 = 1.0 / (GetGamma() - 1.0);
  Thermodynamics* pthermo = pmy_block_->pimpl->pthermo;

  // Force outer-loop vectorization
#pragma omp simd
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      // #pragma omp simd
#pragma novector
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
        for (int n = 1; n <= NVAPOR; ++n) {
          cons(n, k, j, i) = prim(n, k, j, i) * w_d;
          cons(IDN, k, j, i) -= cons(n, k, j, i);
        }

        // momentum
        u_m1 = w_vx * w_d;
        u_m2 = w_vy * w_d;
        u_m3 = w_vz * w_d;

        // total energy
        Real KE = 0.5 * w_d * (w_vx * w_vx + w_vy * w_vy + w_vz * w_vz);
        Real fsig = 1., feps = 1.;
        // vapors
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += prim(n, k, j, i) * (pthermo->GetCvRatio(n) - 1.);
          feps += prim(n, k, j, i) * (1. / pthermo->GetMassRatio(n) - 1.);
        }
        u_e = igm1 * w_p * fsig / feps + KE;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  Thermodynamics* pthermo = pmy_block_->pimpl->pthermo;

  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += prim[n] * (pthermo->GetCvRatio(n) - 1.);
    feps += prim[n] * (1. / pthermo->GetMassRatio(n) - 1.);
  }

  return std::sqrt((1. + (gamma_ - 1) * feps / fsig) * prim[IPR] / prim[IDN]);
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
  Real gm1 = GetGamma() - 1.0;
  Real& w_d = prim(IDN, k, j, i);
  Real& w_p = prim(IPR, k, j, i);

  Real& u_d = cons(IDN, k, j, i);
  Real& u_e = cons(IEN, k, j, i);
  // apply (prim) density floor, without changing momentum or energy
  w_d = (w_d > density_floor_) ? w_d : density_floor_;
  // ensure cons density matches
  u_d = w_d;

  Real e_k = 0.5 * w_d *
             (SQR(prim(IVX, k, j, i)) + SQR(prim(IVY, k, j, i)) +
              SQR(prim(IVZ, k, j, i)));
  // apply pressure floor, correct total energy
  u_e = (w_p > pressure_floor_) ? u_e : ((pressure_floor_ / gm1) + e_k);
  w_p = (w_p > pressure_floor_) ? w_p : pressure_floor_;

  return;
}
