// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <air_parcel.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <configure.hpp>
#include <impl.hpp>

// climath
#include <climath/core.h>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// astro
#include <astro/celestrial_body.hpp>

// harp
#include "harp/radiation.hpp"

Real P0, T0, Omega, grav, gammad, Tmin, radius;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
#ifdef CUBED_SPHERE
  AllocateUserOutputVariables(9);
#else
  AllocateUserOutputVariables(3);
#endif

  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "pres");
#ifdef CUBED_SPHERE
  SetUserOutputVariableName(3, "lat");
  SetUserOutputVariableName(4, "lon");
  SetUserOutputVariableName(5, "vlat");
  SetUserOutputVariableName(6, "vlon");
  SetUserOutputVariableName(7, "w");
  SetUserOutputVariableName(8, "zenith");
#endif
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, P0, k, j, i);
        user_out_var(2, k, j, i) = phydro->w(IPR, k, j, i);
      }

#ifdef CUBED_SPHERE
  auto pexo3 = pimpl->pexo3;
  Real lat, lon;
  Real U, V;
  Direction ray = pimpl->prad->GetRayInput(0);
  Real zenith;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, phydro->w(IVY, k, j, i), phydro->w(IVZ, k, j, i),
                     k, j, i);
        user_out_var(3, k, j, i) = lat;
        user_out_var(4, k, j, i) = lon;
        user_out_var(5, k, j, i) = U;
        user_out_var(6, k, j, i) = V;
        user_out_var(7, k, j, i) = phydro->w(IVX, k, j, i);

        ray = pimpl->planet->ParentZenithAngle(pmy_mesh->time, M_PI / 2. - lat,
                                               lon);
        zenith = std::acos(ray.mu) / M_PI * 180.0;
        user_out_var(8, k, j, i) = zenith;
      }
#endif
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pthermo = Thermodynamics::GetInstance();
  auto prad = pmb->pimpl->prad;

#ifdef CUBED_SPHERE
  auto pexo3 = pmb->pimpl->pexo3;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);

        // coriolis force
        Real f = 2. * Omega * sin(lat);
        Real f2 = 2. * Omega * cos(lat);
        Real U, V;

        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        Real ll_acc_U = f * m3;
        Real ll_acc_V = -f * m2;
        Real acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        du(IM2, k, j, i) += dt * acc2;
        du(IM3, k, j, i) += dt * acc3;
      }
#endif
}

Real TimeStep(MeshBlock *pmb) {
  auto prad = pmb->pimpl->prad;
  Real time = pmb->pmy_mesh->time;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      prad->CalTimeStep(pmb, k, j, pmb->is, pmb->ie);
    }

  return prad->GetTimeStep() * (1. + time / prad->GetRelaxTime());
}

//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // forcing parameters
  Omega = pin->GetReal("problem", "Omega");
  grav = -pin->GetReal("hydro", "grav_acc1");
  T0 = pin->GetReal("problem", "T0");
  P0 = pin->GetReal("problem", "P0");
  gammad = pin->GetReal("hydro", "gamma");
  Tmin = pin->GetReal("problem", "Tmin");
  radius = pin->GetReal("problem", "radius");

  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
  // EnrollUserTimeStepFunction(TimeStep);
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));
  auto pexo3 = pimpl->pexo3;
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);
  // estimate surface temperature and pressure
  // thermodynamic constants
  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real Rd = pthermo->GetRd();
  Real cp = gammad / (gammad - 1.) * Rd;

  // estimate surface temperature and pressure
  Real Ts = T0 - grav / cp * (x1min - radius);
  Real Ps = P0 * pow(Ts / T0, cp / Rd);

  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.SetZero();
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2., "reversible", grav);

      int i = is;
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        //        pthermo->Extrapolate(&air, pcoord->dx1f(i), "isothermal",
        //        grav);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "dry", grav, 1.e-3);
      }
    }

  for (int i = ie - 1; i >= is; --i) {
    std::cout << "i = " << i << " pres = "
              << (phydro->u(IPR, ks, js, i + 1) + phydro->u(IPR, ks, js, i)) /
                     2.
              << " dens = "
              << (phydro->u(IDN, ks, js, i + 1) + phydro->u(IDN, ks, js, i)) /
                     2.
              << std::endl;
  }
}
