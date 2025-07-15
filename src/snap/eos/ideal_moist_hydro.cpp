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
#include <interface/eos.hpp>

// exo3
#include <exo3/exo3.hpp>

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
  auto pmb = pmy_block_;

  auto u = get_all(cons);
  auto w = pmb->pimpl->peos->compute("U->W", {u});

  auto prim_out = get_all(prim);
  prim_out.copy_(w);
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
  auto pmb = pmy_block_;

  auto w = get_all(prim);
  auto u = pmb->pimpl->peos->compute("W->U", {w});

  auto cons_out = get_all(cons);
  cons_out.copy_(u);

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
#pragma omp declare simd simdlen(SIMD_WIDTH) uniform(this)
Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  auto pmb = pmy_block_;
  Real gamma = get_gamma(pmb->pimpl->peos, prim);
  return std::sqrt(gamma * prim[IPR] / prim[IDN]);
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int
// k, int j,
//                                                 =int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface
// states
void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real>& prim, int k,
                                           int j, int i) {}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real>
// &prim,
//           AthenaArray<Real> &cons, FaceField &b, int k, int j, int i) {
// \brief Apply pressure (prim) floor and correct energy (cons) (typically after
// W(U))
void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real>& prim,
                                                    AthenaArray<Real>& cons,
                                                    AthenaArray<Real>& bcc,
                                                    int k, int j, int i) {}
