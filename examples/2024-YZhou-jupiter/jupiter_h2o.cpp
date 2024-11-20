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
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// special includes
// #include <special/giants_enroll_vapor_functions_v1.hpp>

Real grav, P0, T0, Tmin, prad, hrate;
int iH2O;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(4 + NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "rh_H2O");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) = potential_temp(pthermo, w.at(k, j, i), P0);
        // theta_v
        user_out_var(2, k, j, i) =
            user_out_var(1, k, j, i) * pthermo->RovRd(w.at(k, j, i));
        // mse
        user_out_var(3, k, j, i) =
            moist_static_energy(pthermo, w.at(k, j, i), grav * pcoord->x1v(i));
        // relative humidity
        user_out_var(4, k, j, i) = relative_humidity(pthermo, w.at(k, j, i))[1];
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real cv = pthermo->GetCv(w.at(k, j, i));
        if (w(IPR, k, j, i) < prad) {
          du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
                              (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));
        }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");

  Tmin = pin->GetReal("problem", "Tmin");
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate") / 86400.;

  // index
  iH2O = pthermo->SpeciesIndex("H2O");
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: jupiter_crm");

  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

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

  Real Ts = T0 - grav / cp * x1min;
  Real Ps = P0 * pow(Ts / T0, cp / Rd);
  Real yH2O = pin->GetReal("problem", "qH2O.gkg") / 1.E3;

  std::vector<Real> yfrac(IVX, 0.);
  yfrac[iH2O] = yH2O;
  yfrac[0] = 1. - yH2O;

  while (iter++ < max_iter) {
    pthermo->SetMassFractions<Real>(yfrac.data());
    pthermo->EquilibrateTP(Ts, Ps);

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate_inplace(pcoord->dx1f(i), "pseudo", grav);
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
      pthermo->SetMassFractions<Real>(yfrac.data());
      pthermo->EquilibrateTP(Ts, Ps);

      // half a grid to cell center
      pthermo->Extrapolate_inplace(pcoord->dx1f(is) / 2., "reversible", grav);

      int i = is;
      for (; i <= ie; ++i) {
        pthermo->GetPrimitive(w.at(k, j, i));
        if (pthermo->GetTemp() < Tmin) break;
        pthermo->Extrapolate_inplace(pcoord->dx1f(i), "pseudo", grav);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        pthermo->GetPrimitive(w.at(k, j, i));
        pthermo->Extrapolate_inplace(pcoord->dx1f(i), "isothermal", grav);
      }
    }
}
