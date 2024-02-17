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
#include <athena/stride_iterator.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere_utility.hpp>
#include <exo3/gnomonic_equiangle.hpp>

// snap
#include "../thermodynamics/thermodynamics.hpp"
#include "eos_helper.hpp"

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#endif

std::string str_grid_ij(AthenaArray<Real> const& var, int n, int k, int j,
                        int i, int il, int iu, int jl, int ju) {
  std::stringstream msg;
  msg << std::endl << std::setw(8) << " ";
  for (int jj = std::max(j - NGHOST, jl); jj <= std::min(j + NGHOST, ju);
       ++jj) {
    msg << "j = " << jj << ", ";
  }
  msg << std::endl;

  for (int ii = std::max(i - NGHOST, il); ii <= std::min(i + NGHOST, iu);
       ++ii) {
    msg << "i = " << ii << " |";
    for (int jj = std::max(j - NGHOST, jl); jj <= std::min(j + NGHOST, ju);
         ++jj)
      msg << std::setw(8) << var(n, k, jj, ii) << ", ";
    msg << " |" << std::endl;
  }

  return msg.str();
}

namespace cs = CubedSphereUtility;

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
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmy_block_;

  apply_vapor_limiter(&cons, pmy_block_);

  Real gm1 = GetGamma() - 1.0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
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
        for (int n = 0; n <= NVAPOR; ++n) density += cons(n, k, j, i);
        w_d = density;
        Real di = 1. / density;

#ifdef ENABLE_GLOG
        LOG_IF(ERROR, std::isnan(w_d) || (w_d < density_floor_))
            << "cycle = " << pmb->pmy_mesh->ncycle << ", "
            << "rank = " << Globals::my_rank << ", (k,j,i) = "
            << "(" << k << "," << j << "," << i << ")"
            << ", w_d = " << w_d << std::endl;

        LOG_IF(INFO, std::isnan(w_d) || (w_d < density_floor_))
            << str_grid_ij(cons, IDN, k, j, i, il, iu, jl, ju);
#endif

        // mass mixing ratio
        for (int n = 1; n <= NVAPOR; ++n)
          prim(n, k, j, i) = cons(n, k, j, i) * di;

#ifdef ENABLE_GLOG
        LOG_IF(ERROR, std::isnan(w_d) || (w_d < density_floor_))
            << "rank = " << Globals::my_rank << ", (k,j,i) = "
            << "(" << k << "," << j << "," << i << ")"
            << ", w_d = " << w_d << std::endl;

        LOG_IF(INFO, std::isnan(w_d) || (w_d < density_floor_))
            << str_grid_ij(cons, IDN, k, j, i, il, iu, jl, ju);
#endif

        // internal energy
        Real KE, fsig = 1., feps = 1.;
        // vapors
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += prim(n, k, j, i) * (pthermo->GetCvRatioMass(n) - 1.);
          feps += prim(n, k, j, i) * (pthermo->GetInvMuRatio(n) - 1.);
        }

#ifdef CUBED_SPHERE
        Real cos_theta =
            static_cast<GnomonicEquiangle*>(pco)->GetCosineCell(k, j);
#endif

        // Real di = 1.0/u_d;
        w_vx = u_m1 * di;
        w_vy = u_m2 * di;
        w_vz = u_m3 * di;

#ifdef CUBED_SPHERE
        cs::CovariantToContravariant(prim.at(k, j, i), cos_theta);
#endif

        // internal energy
        KE = 0.5 * (u_m1 * w_vx + u_m2 * w_vy + u_m3 * w_vz);
        w_p = gm1 * (u_e - KE) * feps / fsig;

#ifdef ENABLE_GLOG
        LOG_IF(WARNING, (w_p < 0.) || std::isnan(w_p))
            << "cycle = " << pmb->pmy_mesh->ncycle << ", "
            << "rank = " << Globals::my_rank << ", (k,j,i) = "
            << "(" << k << "," << j << "," << i << ")"
            << ", w_p = " << w_p << std::endl;
#endif
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
  Real igm1 = 1.0 / (GetGamma() - 1.0);
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
        for (int n = 1; n <= NVAPOR; ++n) {
          cons(n, k, j, i) = prim(n, k, j, i) * w_d;
          cons(IDN, k, j, i) -= cons(n, k, j, i);
        }

        // momentum
        u_m1 = w_vx * w_d;
        u_m2 = w_vy * w_d;
        u_m3 = w_vz * w_d;

#ifdef CUBED_SPHERE
        cs::ContravariantToCovariant(
            cons.at(k, j, i),
            static_cast<GnomonicEquiangle*>(pco)->GetCosineCell(k, j));
#endif

        // total energy
        Real KE = 0.5 * (u_m1 * w_vx + u_m2 * w_vy + u_m3 * w_vz);
        Real fsig = 1., feps = 1.;
        // vapors
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += prim(n, k, j, i) * (pthermo->GetCvRatioMass(n) - 1.);
          feps += prim(n, k, j, i) * (pthermo->GetInvMuRatio(n) - 1.);
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
  auto pthermo = Thermodynamics::GetInstance();

  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += prim[n] * (pthermo->GetCvRatioMass(n) - 1.);
    feps += prim[n] * (pthermo->GetInvMuRatio(n) - 1.);
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
