//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file hs94.cpp
//  \brief Problem generator for Held-Suarez-94 GCM bench mark.
//
// REFERENCE: I.M Held & M.J Suarez, "A Proposal for the Intercomparison of the
// Dynamical Cores of Atmospheric General Circulation Models"
// C++ headers
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

static Real p0, Omega, Rd, cp, Kf, Ts, Rp, z_iso, sponge_lat, sponge_tau, grav,
    vis, heat_flux;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

// \brief Held-Suarez atmosphere benchmark test. Refernce: Held & Suarez,
// (1994). Forcing parameters are given in the paper.
//! \fn void Damping(MeshBlock *pmb, Real const time, Real const dt,
//  AthenaArray<Real> const& w, AthenaArray<Real> const& bcc, AthenaArray<Real>
//  &u) \brief Pseudo radiative damping of Earth atmosphere for HS94 test.

// The x1x2 is now x2x3
void FindLatlon(Real *lat, Real *lon, Real x2, Real x1) {
  Real dist = sqrt(x1 * x1 + x2 * x2);
  *lat = M_PI / 2. - dist / Rp;
  *lon = asin(x1 / dist);
  if (x2 > 0 && x1 > 0) *lon = M_PI - *lon;
  if (x2 > 0 && x1 < 0) *lon = -M_PI - *lon;
}

void FindXY(Real *x1, Real *x2, Real lat, Real lon) {
  Real dist = Rp * (M_PI / 2. - lat);
  *x1 = dist * sin(lon);
  *x2 = -dist * cos(lon);
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &cons_scalar) {
  // Coriolis force
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real omega1, omega2;
        // Load Latitude
        Real lat, lon;
        FindLatlon(&lat, &lon, pmb->pcoord->x3v(i), pmb->pcoord->x2v(j));

        omega1 = cos(lat) * Omega;
        omega2 = sin(lat) * Omega;  // f
        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * w(IVY, k, j, i);
        Real m3 = w(IDN, k, j, i) * w(IVZ, k, j, i);
        u(IM1, k, j, i) -= -2. * dt * (omega1 * m2);
        u(IM2, k, j, i) -= 2. * dt * (omega1 * m1 - omega2 * m3);
        u(IM3, k, j, i) -= 2. * dt * (omega2 * m2);
      }
    }
  }

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Momentum damping if high in the atmosphere
        if (pmb->pcoord->x1v(i) > z_iso) {
          Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
          Real m2 = w(IDN, k, j, i) * w(IVY, k, j, i);
          Real m3 = w(IDN, k, j, i) * w(IVZ, k, j, i);
          u(IM1, k, j, i) += -dt * Kf * m1;
          u(IM2, k, j, i) += -dt * Kf * m2;
          u(IM3, k, j, i) += -dt * Kf * m3;
        }
      }

  // Treat the boundaries
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      // Energy input 7.7 W/m^2
      Real dz = pmb->pcoord->dx1f(pmb->is);
      u(IEN, k, j, pmb->is) += heat_flux * dt / dz;
      // Energy taken away from the top
      dz = pmb->pcoord->dx1f(pmb->ie);
      u(IEN, k, j, pmb->ie) -= heat_flux * dt / dz;
      // if(u(IEN, k, j, pmb->ie)<0) u(IEN, k, j, pmb->ie) = 0;
    }

  // Sponge layer and viscosity
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real dist = sqrt(pmb->pcoord->x3v(i) * pmb->pcoord->x3v(i) +
                         pmb->pcoord->x2v(j) * pmb->pcoord->x2v(j));
        Real lat_now = 90. - dist / Rp / M_PI * 180.;
        Real s = lat_now / sponge_lat;

        if (s < 1) {
          // Softer sponge layer, linear increasing from 0
          u(IM1, k, j, i) -= dt * w(IDN, k, j, i) * w(IVX, k, j, i) /
                             sponge_tau * (sponge_lat - lat_now) / 5;
          u(IM2, k, j, i) -= dt * w(IDN, k, j, i) * w(IVY, k, j, i) /
                             sponge_tau * (sponge_lat - lat_now) / 5;
          u(IM3, k, j, i) -= dt * w(IDN, k, j, i) * w(IVZ, k, j, i) /
                             sponge_tau * (sponge_lat - lat_now) / 5;
        } else {  // viscosity
          Real dz = pmb->pcoord->dx1f(i);
          Real dx = pmb->pcoord->dx2f(j);
          Real dy = pmb->pcoord->dx3f(k);
          Real dir = IVX;
          Real lap_1 = (w(dir, k, j, i + 1) - 2 * w(dir, k, j, i) +
                        w(dir, k, j, i - 1)) /
                           (dz * dz) +
                       (w(dir, k, j + 1, i) - 2 * w(dir, k, j, i) +
                        w(dir, k, j - 1, i)) /
                           (dx * dx) +
                       (w(dir, k + 1, j, i) - 2 * w(dir, k, j, i) +
                        w(dir, k - 1, j, i)) /
                           (dy * dy);
          dir = IVY;
          Real lap_2 = (w(dir, k, j, i + 1) - 2 * w(dir, k, j, i) +
                        w(dir, k, j, i - 1)) /
                           (dz * dz) +
                       (w(dir, k, j + 1, i) - 2 * w(dir, k, j, i) +
                        w(dir, k, j - 1, i)) /
                           (dx * dx) +
                       (w(dir, k + 1, j, i) - 2 * w(dir, k, j, i) +
                        w(dir, k - 1, j, i)) /
                           (dy * dy);
          dir = IVZ;
          Real lap_3 = (w(dir, k, j, i + 1) - 2 * w(dir, k, j, i) +
                        w(dir, k, j, i - 1)) /
                           (dz * dz) +
                       (w(dir, k, j + 1, i) - 2 * w(dir, k, j, i) +
                        w(dir, k, j - 1, i)) /
                           (dx * dx) +
                       (w(dir, k + 1, j, i) - 2 * w(dir, k, j, i) +
                        w(dir, k - 1, j, i)) /
                           (dy * dy);
          u(IM1, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_1;
          u(IM2, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_2;
          u(IM3, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_3;
        }
      }
    }
  }
}
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real day_to_s = 8.64E4;
  // forcing parameters
  Omega = pin->GetReal("problem", "Omega");
  // thermodynamic parameters
  Real gamma = pin->GetReal("hydro", "gamma");
  grav = pin->GetReal("hydro", "grav_acc1");
  Ts = pin->GetReal("problem", "Ts");
  p0 = pin->GetReal("problem", "p0");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  // damping parameters
  Kf = pin->GetReal("problem", "Kf");
  Kf /= day_to_s;
  z_iso = pin->GetReal("problem", "z_iso");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  vis = pin->GetReal("problem", "vis");
  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
}
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
// Problem generator kept the same for init conditions
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  grav = pin->GetReal("hydro", "grav_acc1");
  Real gamma = pin->GetReal("hydro", "gamma");
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  Rp = pin->GetReal("problem", "Rp");
  z_iso = pin->GetReal("problem", "z_iso");
  sponge_lat = pin->GetReal("problem", "sponge_lat");
  heat_flux = pin->GetReal("problem", "heat_flux");
  // construct an adiabatic atmosphere
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[IPR] = p0;
      air.w[IDN] = Ts;

      int i = is;
      // for (; i <= ie; ++i) {
      //   if (pcoord->x1v(i) - Rp > z_iso) break;

      //   pimpl->DistributeToConserved(air, k, j, i);
      //   pthermo->Extrapolate(&air, pcoord->dx1f(i),
      //                        Thermodynamics::Method::DryAdiabat, grav,
      //                        0.001);
      //   // add noise
      //   // air.w[IVX] = 0.01 * distribution(generator);
      // }

      // construct isothermal atmosphere
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::Isothermal, grav, 0.001);
        // add noise
        air.w[IVY] = 0.00001 * distribution(generator);
        air.w[IVZ] = 0.00001 * distribution(generator);
      }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

// \brif Output distributions of temperature and potential temperature.
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real prim[NHYDRO];
        for (int n = 0; n < NHYDRO; ++n) prim[n] = phydro->w(n, j, i);
        Real temp = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / Rd;
        user_out_var(0, k, j, i) = temp;
        user_out_var(1, k, j, i) =
            temp * pow(p0 / phydro->w(IPR, k, j, i), Rd / cp);
      }
}
