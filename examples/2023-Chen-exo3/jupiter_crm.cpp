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
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>


#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

Real grav, P0, T0, Tmin, prad, hrate, Omega, Rp;
int iH2O, iNH3;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

namespace cs = CubedSphereUtility;

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();

  auto pexo3 = pmb->pimpl->pexo3;
  Real om_earth = Omega;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }


  /* Sponge Layer
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real pres = w(IPR, k, j, i);
        if (pres < piso) {  // sponge layer at top
          Real tau = sponge_tau * pow(pres / piso, 2);
          u(IVX, k, j, i) = u(IVX, k, j, i) / (1 + dt / tau);
          u(IVY, k, j, i) = u(IVY, k, j, i) / (1 + dt / tau);
          u(IVZ, k, j, i) = u(IVZ, k, j, i) / (1 + dt / tau);
        }
      }
    }
  }*/
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // coriolis force
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
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(4 + NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3 + n, name.c_str());
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, P0, k, j, i);
        // theta_v
        user_out_var(2, k, j, i) =
            user_out_var(1, k, j, i) * pthermo->RovRd(this, k, j, i);
        // mse
        user_out_var(3, k, j, i) =
            pthermo->MoistStaticEnergy(this, grav * (pcoord->x1v(i)-Rp), k, j, i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3 + n, k, j, i) =
              pthermo->RelativeHumidity(this, n, k, j, i);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Omega = pin->GetReal("problem", "Omega");
  Rp = pin->GetReal("problem", "Rp");

  Tmin = pin->GetReal("problem", "Tmin");
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate") / 86400.;

  // index
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVaporId("H2O");
  iNH3 = pindex->GetVaporId("NH3");
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: jupiter_crm");

  auto pthermo = Thermodynamics::GetInstance();

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  // request temperature and pressure
  app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;

  // set up an adiabatic atmosphere
  int max_iter = 400, iter = 0;
  Real Ttol = pin->GetOrAddReal("problem", "init_Ttol", 0.01);

  AirParcel air(AirParcel::Type::MoleFrac);

  // estimate surface temperature and pressure
  Real Ts = T0 - grav / cp * (x1min - Rp);
  Real Ps = P0 * pow(Ts / T0, cp / Rd);
  Real xH2O = pin->GetReal("problem", "qH2O.ppmv") / 1.E6;
  Real xNH3 = pin->GetReal("problem", "qNH3.ppmv") / 1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    air.w[iH2O] = xH2O;
    air.w[iNH3] = xNH3;
    air.w[IPR] = Ps;
    air.w[IDN] = Ts;

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav);
      if (air.w[IPR] < P0) break;
    }

    // make up for the difference
    Ts += T0 - air.w[IDN];
    if (std::abs(T0 - air.w[IDN]) < Ttol) break;

    app->Log("Iteration #", iter);
    app->Log("T", air.w[IDN]);
  }

  if (iter > max_iter) {
    throw RuntimeError("ProblemGenerator", "maximum iteration reached");
  }

  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.SetZero();
      air.w[iH2O] = xH2O;
      air.w[iNH3] = xNH3;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2., "reversible", grav);

      int i = is;
      for (; i <= ie; ++i) {
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav, 1.e-5);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "isothermal", grav);
      }
    }
}
