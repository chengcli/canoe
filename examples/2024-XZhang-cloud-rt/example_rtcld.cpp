// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

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
#include "svp_functions.hpp"

// astro
#include <astro/celestrial_body.hpp>

// harp
#include "harp/radiation.hpp"

#define _sqr(x) ((x) * (x))
#define _qur(x) ((x) * (x) * (x) * (x))

using namespace std;

Real grav, P0, T0, Tmin, Omega, radius;
int iH2O, iNH3;

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 1.0);

namespace cs = CubedSphereUtility;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
#ifdef CUBED_SPHERE
  AllocateUserOutputVariables(11 + NVAPOR);
#else
  AllocateUserOutputVariables(5 + NVAPOR);
#endif

  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "pres");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(4 + n, name.c_str());
  }

#ifdef CUBED_SPHERE
  SetUserOutputVariableName(5 + NVAPOR, "lat");
  SetUserOutputVariableName(6 + NVAPOR, "lon");
  SetUserOutputVariableName(7 + NVAPOR, "vlat");
  SetUserOutputVariableName(8 + NVAPOR, "vlon");
  SetUserOutputVariableName(9 + NVAPOR, "w");
  SetUserOutputVariableName(10 + NVAPOR, "zenith");
#endif
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
        user_out_var(3, k, j, i) = pthermo->MoistStaticEnergy(
            this, grav * (pcoord->x1v(i) - radius), k, j, i);
        user_out_var(4, k, j, i) = phydro->w(IPR, k, j, i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(4 + n, k, j, i) =
              pthermo->RelativeHumidity(this, n, k, j, i);
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
        user_out_var(5 + NVAPOR, k, j, i) = lat;
        user_out_var(6 + NVAPOR, k, j, i) = lon;
        user_out_var(7 + NVAPOR, k, j, i) = U;
        user_out_var(8 + NVAPOR, k, j, i) = V;
        user_out_var(9 + NVAPOR, k, j, i) = phydro->w(IVX, k, j, i);

        ray = pimpl->planet->ParentZenithAngle(pmy_mesh->time, M_PI / 2. - lat,
                                               lon);
        zenith = std::acos(ray.mu) / M_PI * 180.0;
        user_out_var(10 + NVAPOR, k, j, i) = zenith;
      }
#endif
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();

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

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Omega = pin->GetReal("problem", "Omega");
  radius = pin->GetReal("problem", "radius");
  Tmin = pin->GetReal("problem", "Tmin");

  // index
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVaporId("H2O");
  iNH3 = pindex->GetVaporId("NH3");
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));
  auto pthermo = Thermodynamics::GetInstance();
  AirParcel air(AirParcel::Type::MoleFrac);

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  // thermodynamic constants
  Real gammad = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gammad / (gammad - 1.) * Rd;

  // set up an adiabatic atmosphere
  int max_iter = 400, iter = 0;
  Real Ttol = pin->GetOrAddReal("problem", "init_Ttol", 0.01);

  // estimate surface temperature and pressure
  Real Ts = T0 - grav / cp * (x1min - radius);
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
