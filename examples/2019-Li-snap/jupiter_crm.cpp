/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2023
 * Contact: chengcli@umich.edu
 * Reference: Test Jupiter CRM
 * -------------------------------------------------------------------------------------
 */

// athena
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>
#include <index_map.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>

Real grav, P0, T0, Tmin, prad, hrate;
int iH2O, iNH3;

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

void FindLatlon(Real *lat, Real *lon, Real x2, Real x1) {
  Real Rp = 7E7;
  Real dist = sqrt(x1 * x1 + x2 * x2);
  *lat = M_PI / 2. - dist / Rp;
  *lon = asin(x1 / dist);
  if (x2 > 0 && x1 > 0) *lon = M_PI - *lon;
  if (x2 > 0 && x1 < 0) *lon = -M_PI - *lon;
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
            pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3 + n, k, j, i) =
              pthermo->RelativeHumidity(this, n, k, j, i);
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();

  // Cooling and Heating
  /*for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);

        Real cv = pthermo->GetCvMass(air, 0);

        if (w(IPR, k, j, i) < prad) {
          du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
                              (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));
        }

        // if (air.w[IDN] < Tmin) {
        //   u(IEN,k,j,i) += w(IDN,k,j,i)*cv*(Tmin - temp)/sponge_tau*dt;
        // }
      }
  */
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real dz = pmb->pcoord->dx1f(pmb->is);
      du(IEN, k, j, is) +=
          hrate * dt / dz * (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));
      ;
      // Energy taken away from the top
      dz = pmb->pcoord->dx1f(pmb->ie);
      for (int i = is; i <= ie; ++i) {
        if (w(IPR, k, j, i) < prad) {  // relax to Tmin profile
          auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);
          Real cv = pthermo->GetCvMass(air, 0);
          Real Rd = pthermo->GetRd();
          Real Tnow = w(IPR, k, j, i) / Rd / w(IDN, k, j, i);
          // du(IEN, k, j, i) += (Tmin - Tnow) * cv * w(IDN, k, j, i) * dt /
          // 1E5;
        }
      }
      // du(IEN, k, j, ie) -= hrate * dt / dz * (1. + 1.E-4 * sin(2. * M_PI *
      // rand() / RAND_MAX));;
    }

  // Coriolis force
  Real Omega = 1.74E-4;
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
        du(IM1, k, j, i) -= -2. * dt * (omega1 * m2);
        du(IM2, k, j, i) -= 2. * dt * (omega1 * m1 - omega2 * m3);
        du(IM3, k, j, i) -= 2. * dt * (omega2 * m2);
      }
    }
  }
  // Sponge layer and viscosity
  Real sponge_lat = 75.0;
  Real sponge_tau = 1E4;
  Real vis = 0.;
  Real Rp = 7E7;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real dist = sqrt(pmb->pcoord->x3v(i) * pmb->pcoord->x3v(i) +
                         pmb->pcoord->x2v(j) * pmb->pcoord->x2v(j));
        Real lat_now = 90. - dist / Rp / M_PI * 180.;
        Real s = lat_now / sponge_lat;

        if (s < 1) {
          // Softer sponge layer, linear increasing from 0
          du(IM1, k, j, i) -= dt * w(IDN, k, j, i) * w(IVX, k, j, i) /
                              sponge_tau * (sponge_lat - lat_now) / 5;
          du(IM2, k, j, i) -= dt * w(IDN, k, j, i) * w(IVY, k, j, i) /
                              sponge_tau * (sponge_lat - lat_now) / 5;
          du(IM3, k, j, i) -= dt * w(IDN, k, j, i) * w(IVZ, k, j, i) /
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
          du(IM1, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_1;
          du(IM2, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_2;
          du(IM3, k, j, i) += vis * dt * w(IDN, k, j, i) * lap_3;
        }
      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");

  // Real Omega = pin->GetReal("problem", "Omega");
  // Real Rp = pin->GetReal("problem", "Rp");

  Tmin = pin->GetReal("problem", "Tmin");
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate") / 86400.;

  // index
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVaporId("H2O");
  // iNH3 = pindex->GetVaporId("NH3");
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
  Real Ts = T0 - grav / cp * x1min;
  Real Ps = P0 * pow(Ts / T0, cp / Rd);
  Real xH2O = pin->GetReal("problem", "qH2O.ppmv") / 1.E6;
  Real xNH3 = pin->GetReal("problem", "qNH3.ppmv") / 1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    air.w[iH2O] = xH2O;
    // air.w[iNH3] = xNH3;
    air.w[IPR] = Ps;
    air.w[IDN] = Ts;

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate(&air, pcoord->dx1f(i),
                           Thermodynamics::Method::PseudoAdiabat, grav);
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
      // air.w[iNH3] = xNH3;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2.,
                           Thermodynamics::Method::ReversibleAdiabat, grav);

      int i = is;
      for (; i <= ie; ++i) {
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::PseudoAdiabat, grav,
                             1.e-5);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::Isothermal, grav);
      }
    }
}
