//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file hjupiter.cpp
//  \brief Problem generator for Hot Juiter GCM benchmark test.
//
// REFERENCE: I.M Held & M.J Suarez, "A Proposal for the Intercomparison of the
// Dynamical Cores of Atmospheric General Circulation Models"

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

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

static Real p0, Omega, Rd, cp, Ts, dT_stra, Gamma_trop, Kt, Rp, beta_trop,
    dT_e2p, sponge_tau, sponge_width, grav, z_stra, sigma_stra, T_eq, T_vert,
    x1min;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

// \brief Held-Suarez atmosphere benchmark test. Refernce: Held & Suarez,
// (1994). Forcing parameters are given in the paper.

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real om_earth = Omega;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        // coriolis force
        Real f = 2. * om_earth * sin(lat);
        Real f2 = 2. * om_earth * cos(lat);
        Real U, V;

        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);

        Real m1 = w(IDN, k, j, i) * w(IVX, k, j, i);
        Real m2 = w(IDN, k, j, i) * U;
        Real m3 = w(IDN, k, j, i) * V;

        // du(IM1, k, j, i) += dt * f * m2;
        Real ll_acc_U = f * m3;  //- f2 * m1;
        Real ll_acc_V = -f * m2;
        Real acc1, acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        du(IM2, k, j, i) += dt * acc2;
        du(IM3, k, j, i) += dt * acc3;
      }

  Real kappa;  // Rd/Cp
  kappa = Rd / cp;
  // Real piso = p0*pow(iso_temp/Ts, cp/Rd);

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Load height, latitude, and longitude
        Real z = pmb->pcoord->x1v(i) - Rp;
        Real theta, phi;
        pexo3->GetLatLon(&theta, &phi, k, j, i);
        // Set-up sigma coordinate.
        // Real T_now = w(IPR,k,j,i)/w(IDN,k,j,i)/Rd;
        // Real p0_now = w(IPR, k, j,
        // 0)*exp((pmb->pcoord->x1v(0)-Rp)*(-grav)/(Rd*T_now));
        Real p0_now = w(IPR, k, j, 0);
        Real sigma = w(IPR, k, j, i) / p0_now;

        if (z <= z_stra) {
          beta_trop =
              sin(M_PI * (sigma - sigma_stra) / (2. * (1. - sigma_stra)));
        } else {
          beta_trop = 0.;
        }

        // Temperature (energy) damping coefficient
        // Temperature difference, T - Teq
        if (z <= z_stra) {
          T_vert = Ts - Gamma_trop * ((z_stra + z) / 2.) +
                   sqrt(_sqr(Gamma_trop * (z - z_stra) / 2.) + _sqr(dT_stra));
        } else {
          T_vert = Ts - Gamma_trop * z_stra + dT_stra;
        }
        Real temp = w(IPR, k, j, i) / w(IDN, k, j, i) / Rd;
        // Diabatic Cooling
        T_eq = T_vert + beta_trop * dT_e2p * cos(theta) * cos(phi);
        du(IEN, k, j, i) +=
            -dt * (cp - Rd) * w(IDN, k, j, i) * (temp - T_eq) / Kt;
      }

  // Sponge Layer
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real z = pmb->pcoord->x1v(i) - Rp;
        if (z > 2.5E6) {  // sponge layer at top and bottom
          Real tau = sponge_tau;
          du(IVX, k, j, i) -= w(IVX, k, j, i) * (dt / tau) * w(IDN, k, j, i);
          du(IVY, k, j, i) -= w(IVY, k, j, i) * (dt / tau) * w(IDN, k, j, i);
          du(IVZ, k, j, i) -= w(IVZ, k, j, i) * (dt / tau) * w(IDN, k, j, i);
        }
      }
    }
  }
}

Real AngularMomentum(MeshBlock *pmb, int iout) {
  auto pexo3 = pmb->pimpl->pexo3;
  Real AMz = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks,
      ke = pmb->ke;

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        Real x1l = pmb->pcoord->x1f(i);
        Real x1u = pmb->pcoord->x1f(i + 1);
        Real U, V;
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, pmb->phydro->w(IVY, k, j, i),
                     pmb->phydro->w(IVZ, k, j, i), k, j, i);

        Real xt = tan(pmb->pcoord->x2v(j));
        Real yt = tan(pmb->pcoord->x3v(k));
        Real sin_theta =
            sqrt((1.0 + xt * xt + yt * yt) / (1.0 + xt * xt) / (1.0 + yt * yt));

        Real x1 = tan(pmb->pcoord->x2f(j));
        Real x2 = tan(pmb->pcoord->x2f(j + 1));
        Real y = tan(pmb->pcoord->x3v(k));
        Real delta1 = sqrt(1.0 + x1 * x1 + y * y);
        Real delta2 = sqrt(1.0 + x2 * x2 + y * y);
        Real dx2_ang = acos(1 / (delta1 * delta2) * (1 + x1 * x2 + y * y));

        Real x = tan(pmb->pcoord->x2v(j));
        Real y1 = tan(pmb->pcoord->x3f(k));
        Real y2 = tan(pmb->pcoord->x3f(k + 1));
        delta1 = sqrt(1.0 + x * x + y1 * y1);
        delta2 = sqrt(1.0 + x * x + y2 * y2);
        Real dx3_ang = acos(1 / (delta1 * delta2) * (1 + x * x + y1 * y2));

        Real vol = pmb->pcoord->dx1f(i) * dx2_ang * dx3_ang * sin_theta;

        // Originally here used cos(lat), which is x2v-pi, strange
        AMz += pmb->phydro->w(IDN, k, j, i) * vol *
               sqrt((_sqr(x1l) + _sqr(x1u)) / 2.) * cos(lat) *
               (Omega * sqrt(0.5 * (_sqr(x1l) + _sqr(x1u))) * cos(lat) + U);
      }
    }
  }

  return AMz;
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
  dT_e2p = pin->GetReal("thermodynamics", "dT_e2p");
  dT_stra = pin->GetReal("thermodynamics", "dT_stra");
  z_stra = pin->GetReal("thermodynamics", "z_stra");
  Gamma_trop = pin->GetReal("thermodynamics", "Gamma_trop");
  Rd = pin->GetReal("thermodynamics", "Rd");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  x1min = pin->GetReal("mesh", "x1min");
  cp = gamma / (gamma - 1.) * Rd;
  // damping parameters
  Kt = pin->GetReal("problem", "Kt");
  sigma_stra = pin->GetReal("problem", "sigma_stra");
  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
  AllocateUserHistoryOutput(1);
  EnrollUserHistoryOutput(0, AngularMomentum, "z-angular-mom");
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  Real grav = -pin->GetReal("hydro", "grav_acc1");
  Real gamma = pin->GetReal("hydro", "gamma");
  p0 = pin->GetReal("problem", "p0");
  Ts = pin->GetReal("problem", "Ts");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma / (gamma - 1.) * Rd;
  Rp = pin->GetReal("problem", "Rp");

  // construct an isothermal atmosphere
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[IPR] = p0;
      air.w[IDN] = Ts;

      int i = is;
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::Isothermal, grav, 0.001);
        // add noise
        air.w[IVY] = 10. * distribution(generator);
        air.w[IVZ] = 10. * distribution(generator);
      }
    }

  // Set up perturbation as the initial condition to break the symmetry of the
  // original initial condition.
  // Real k3 = 5.;

  // for (int k = ks; k <= ke; ++k) {
  //   for (int j = js; j <= je; ++j) {
  //     for (int i = is; i <= ie; ++i) {
  //       Real theta, phi;
  //       pexo3->GetLatLon(&theta, &phi, k, j, i);
  //       Real z = pcoord->x1v(i)-Rp;
  //       // Set-up isothermal atmosphere
  //       Real temp = Ts;
  //       phydro->w(IPR,k,j,i) = p0*exp(-z*grav/(temp*Rd));
  //       phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)/
  //            (Rd* ( temp + 20.*(distribution(generator) - 0.5) * (1. +
  //            cos(k3*phi)) ) );
  //       //phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)/(Rd*temp);
  //       phydro->w(IM1,k,j,i) = 0.;
  //       phydro->w(IM2,k,j,i) = 0.;
  //       phydro->w(IM3,k,j,i) = 0.;
  //     }
  //   }
  // }

  // transfer to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  // peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is,
  // ie, js, je, ks, ke);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "lat");
  SetUserOutputVariableName(3, "lon");
  SetUserOutputVariableName(4, "vlat");
  SetUserOutputVariableName(5, "vlon");
  SetUserOutputVariableName(6, "Teq");
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
        Real lat, lon;
        Real U, V;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        pexo3->GetUV(&U, &V, phydro->w(IVY, k, j, i), phydro->w(IVZ, k, j, i),
                     k, j, i);
        user_out_var(2, k, j, i) = lat;
        user_out_var(3, k, j, i) = lon;
        user_out_var(4, k, j, i) = U;
        user_out_var(5, k, j, i) = V;

        Real kappa;  // Rd/Cp
        kappa = Rd / cp;
        Real z = pcoord->x1v(i) - Rp;
        Real sigma = phydro->w(IPR, k, j, i) / p0;
        if (z <= z_stra) {
          beta_trop =
              sin(M_PI * (sigma - sigma_stra) / (2. * (1. - sigma_stra)));
        } else {
          beta_trop = 0.;
        }
        Real scaled_z = phydro->w(IPR, k, j, i) / p0;
        Real Teq = T_vert + beta_trop * dT_e2p * cos(lat) * cos(lon);
        user_out_var(6, k, j, i) = Teq;
      }
}
