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

// torch
#include <torch/torch.h>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere_utility.hpp>
#include <exo3/gnomonic_equiangle.hpp>

// snap
#include <snap/athena_arrays.hpp>
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

#include "eos.hpp"
#include "eos_helper.hpp"

// checks
#include <checks.hpp>

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
  cons.toDevice(torch::kCPU);
  auto gammad = torch::full_like(cons.tensor()[0], gamma_);

  prim.toDevice(torch::kCPU);
  prim.tensor() = eos_cons2prim_hydro_ideal(cons.tensor(), gammad);
  prim.fromDevice();
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
  auto primc = const_cast<AthenaArray<Real>&>(prim);

  primc.toDevice(torch::kCPU);
  auto gammad = torch::full_like(primc.tensor()[0], gamma_);

  cons.toDevice(torch::kCPU);
  cons.tensor() = eos_prim2cons_hydro_ideal(primc.tensor(), gammad);
  cons.fromDevice();
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
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
