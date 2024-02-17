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
#include <configure.hpp>
#include <impl.hpp>

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

static Real p0, Omega, Rd, cp, sigmab, Kf, Ts, dT, dtheta, Ka, Ks, Rp, scaled_z,
    z_iso, sponge_tau, sponge_width, grav;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

// \brief Held-Suarez atmosphere benchmark test. Refernce: Held & Suarez,
// (1994). Forcing parameters are given in the paper.
//! \fn void Damping(MeshBlock *pmb, Real const time, Real const dt,
//  AthenaArray<Real> const& w, AthenaArray<Real> const& bcc, AthenaArray<Real>
//  &u) \brief Pseudo radiative damping of Earth atmosphere for HS94 test.
void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &cons_scalar) {
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real omega1, omega2;
        // Load Latitude
        Real theta = pmb->pcoord->x2v(j);
        theta = theta - M_PI / 2.;

        omega1 = cos(theta) * Omega;
        omega2 = sin(theta) * Omega;
        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * w(IVY, k, j, i);
        Real m3 = w(IDN, k, j, i) * w(IVZ, k, j, i);
        u(IM1, k, j, i) -= -2. * dt * (omega1 * m3);
        u(IM2, k, j, i) -= 2. * dt * (omega2 * m3);
        u(IM3, k, j, i) -= 2. * dt * (omega1 * m1 - omega2 * m2);
      }
    }
  }

  Real kappa;  // Rd/Cp
  kappa = Rd / cp;
  Real iso_temp = Ts + grav * z_iso / cp;
  // Real piso = p0*pow(iso_temp/Ts, cp/Rd);
  // Newtonian cooling and Rayleigh drag
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Load latitude
        Real theta = pmb->pcoord->x2v(j);
        // Momentum damping coefficient, Kv
        Real scaled_z = w(IPR, k, j, i) / w(IPR, k, j, pmb->is);
        Real sigma = w(IPR, k, j, i) / p0;
        Real sigma_p = (sigma - sigmab) / (1. - sigmab);
        Real Kv = (sigma_p < 0.0) ? 0.0 : sigma_p * Kf;
        // Temperature (energy) damping coefficient
        // Temperature difference, T - Teq
        Real Teq_p = Ts - dT * _sqr(cos(theta)) -
                     dtheta * log(scaled_z) * _sqr(sin(theta));
        Teq_p *= pow(scaled_z, kappa);
        Real Teq = (Teq_p > 200.) ? Teq_p : 200.;
        Real temp =
            pmb->phydro->w(IPR, k, j, i) / pmb->phydro->w(IDN, k, j, i) / Rd;

        // Temperature damping coefficient, Kt
        sigma_p = (sigma_p < 0.0) ? 0.0 : sigma_p * _qur(sin(theta));
        Real Kt = Ka + (Ks - Ka) * sigma_p;
        // Momentum and energy damping
        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * w(IVY, k, j, i);
        Real m3 = w(IDN, k, j, i) * w(IVZ, k, j, i);
        u(IM1, k, j, i) += -dt * Kv * m1;
        u(IM2, k, j, i) += -dt * Kv * m2;
        u(IM3, k, j, i) += -dt * Kv * m3;
        u(IEN, k, j, i) +=
            -dt * (cp - Rd) * w(IDN, k, j, i) * Kt * (temp - Teq);
      }
  /* Sponge layer
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real pres = w(IPR,k,j,i);
        if (pres < piso) {  // sponge layer at top
            Real tau = sponge_tau*pow(pres/piso,2);
            u(IVX,k,j,i) = u(IVX,k,j,i)/(1 + dt/tau);
            u(IVY,k,j,i) = u(IVY,k,j,i)/(1 + dt/tau);
            u(IVZ,k,j,i) = u(IVZ,k,j,i)/(1 + dt/tau);
        }
      }
    }
  }*/
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
  dtheta = pin->GetReal("thermodynamics", "dtheta");
  dT = pin->GetReal("thermodynamics", "dT");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  // damping parameters
  Kf = pin->GetReal("problem", "Kf");
  Kf /= day_to_s;
  Ka = pin->GetReal("problem", "Ka");
  Ka /= day_to_s;
  Ks = pin->GetReal("problem", "Ks");
  Ks /= day_to_s;
  sigmab = pin->GetReal("problem", "sigmab");
  z_iso = pin->GetReal("problem", "z_iso");
  // sponge_tau = pin->GetReal("problem", "sponge_tau");
  sponge_width = pin->GetReal("hydro", "osponge1_width");
  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
}
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  grav = pin->GetReal("hydro", "grav_acc1");
  Real gamma = pin->GetReal("hydro", "gamma");
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  Rp = pin->GetReal("problem", "Rp");
  z_iso = pin->GetReal("problem", "z_iso");
  // construct an adiabatic atmosphere
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[IPR] = p0;
      air.w[IDN] = Ts;

      int i = is;
      for (; i <= ie; ++i) {
        if (pcoord->x1v(i) - Rp > z_iso) break;

        pimpl->DistributeToConserved(air, k, j, i);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::DryAdiabat, grav);
        // add noise
        air.w[IVX] = 0.01 * distribution(generator);
      }

      // construct isothermal atmosphere
      for (; i <= ie; ++i) {
        pimpl->DistributeToConserved(air, k, j, i);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::Isothermal, grav);
      }
    }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}
// \brif Output distributions of temperature and potential temperature.
void MeshBlock::UserWorkInLoop() {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real prim[NHYDRO];
        for (int n = 0; n < NHYDRO; ++n) prim[n] = phydro->w(n, j, i);
        Real temp = phydro->w(IPR, k, j, i) / phydro->w(IDN, k, j, i) / Rd;
        user_out_var(0, k, j, i) = temp;
        user_out_var(1, k, j, i) =
            temp * pow(p0 / phydro->w(IPR, k, j, i), Rd / cp);
      }
    }
  }
}
