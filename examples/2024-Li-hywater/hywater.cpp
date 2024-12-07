/* -------------------------------------------------------------------------------------
 * HYWATER: A Hydrogen-water Cloud Resolving Model
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2024
 * Contact: chengcli@umich.edu
 * Reference: TBD
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

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

Real grav, P0, T0;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  AllocateUserOutputVariables(8 + NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "tempv");
  SetUserOutputVariableName(2, "enthalpy");
  SetUserOutputVariableName(3, "entropy");
  SetUserOutputVariableName(4, "intEng");

  SetUserOutputVariableName(5, "theta");
  SetUserOutputVariableName(6, "thetav");
  SetUserOutputVariableName(7, "mse");

  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + pthermo->SpeciesName(n);
    SetUserOutputVariableName(7 + n, name.c_str());
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) =
            user_out_var(0, k, j, i) * pthermo->RovRd(w.at(k, j, i));
        user_out_var(2, k, j, i) = pthermo->GetEnthalpy(w.at(k, j, i));
        user_out_var(3, k, j, i) = pthermo->GetEntropy(w.at(k, j, i));
        user_out_var(4, k, j, i) = pthermo->GetInternalEnergy(w.at(k, j, i));

        // theta
        user_out_var(5, k, j, i) = potential_temp(pthermo, w.at(k, j, i), P0);
        // theta_v
        user_out_var(6, k, j, i) =
            user_out_var(5, k, j, i) * pthermo->RovRd(w.at(k, j, i));
        // mse
        user_out_var(7, k, j, i) =
            moist_static_energy(pthermo, w.at(k, j, i), grav * pcoord->x1v(i));
        auto rh = relative_humidity(pthermo, w.at(k, j, i));
        for (int n = 1; n <= NVAPOR; ++n) user_out_var(7 + n, k, j, i) = rh[n];
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");

  auto pthermo = Thermodynamics::GetInstance();
  auto kinetics = get_kinetics_object(pthermo);
  kinetics->setBoiling(true);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: jcrm");

  auto pthermo = Thermodynamics::GetInstance();

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  // request temperature and pressure
  app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real Rd = pthermo->GetRd();
  Real cp = pthermo->GetGammad() / (pthermo->GetGammad() - 1.) * Rd;

  // set up atmosphere
  Real Ttol = pin->GetOrAddReal("problem", "Ttol_abs", 0.01);
  Real Ptol = pin->GetOrAddReal("problem", "Ptol_rel", 1.e-4);
  Real adTdz = pin->GetOrAddReal("problem", "adTdz", 0.);
  Real Tmin = pin->GetOrAddReal("problem", "Tmin", 110.);

  // estimate surface temperature and pressure
  Real Ts = T0 - grav / cp * x1min;
  Real Ps = P0 * pow(Ts / T0, cp / Rd);

  // read deep composition
  std::vector<Real> yfrac(IVX, 0.);
  Real qdry = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "q" + pthermo->SpeciesName(n);
    yfrac[n] = pin->GetReal("problem", name);
    qdry -= yfrac[n];
  }
  yfrac[0] = qdry;

  // construct 1d atmosphere from bottom up
  int nx1 = pmy_mesh->mesh_size.nx1;
  Real dz = (x1max - x1min) / (nx1 - 1);

  AthenaArray<Real> w1, z1;
  w1.NewAthenaArray(NHYDRO, nx1);

  z1.NewAthenaArray(nx1);
  z1(0) = x1min;
  for (int i = 1; i < nx1; ++i) z1(i) = z1(i - 1) + dz;

  int max_iter = 200, iter = 0;
  Real prim[NHYDRO], prim1[NHYDRO];
  while (iter++ < max_iter) {
    pthermo->SetMassFractions<Real>(yfrac.data());
    pthermo->EquilibrateTP(Ts, Ps);

    // stop at just above Z = 0
    int i = 0;
    for (; i < nx1; ++i) {
      pthermo->GetPrimitive<Real>(prim1);
      pthermo->Extrapolate_inplace(dz, "pseudo", grav);
      if (z1(i + 1) > 0.) break;
    }
    pthermo->GetPrimitive<Real>(prim);

    // linear interpolate to Z = 0
    for (int n = 0; n < NHYDRO; ++n)
      prim[n] = prim[n] + (prim1[n] - prim[n]) * z1(i + 1) / dz;

    // make up for the difference
    Real t0 = pthermo->GetTemp(prim);
    Real p0 = prim[IPR];
    Ts += T0 - t0;
    Ps *= P0 / p0;
    if ((fabs(T0 - t0) < Ttol) && (fabs(P0 / p0 - 1.) < Ptol)) break;

    app->Log("Iteration #", iter);
    app->Log("T", t0);
    app->Log("p", p0);
  }

  if (iter > max_iter) {
    throw RuntimeError("ProblemGenerator", "maximum iteration reached");
  }

  std::cout << "########## Finished ##########" << std::endl;
  std::cout << "Ts = " << Ts << std::endl;
  std::cout << "Ps = " << Ps << std::endl;

  // construct 1d atmosphere
  pthermo->SetMassFractions<Real>(yfrac.data());
  pthermo->EquilibrateTP(Ts, Ps);

  // half a grid to cell center
  pthermo->Extrapolate_inplace(dz / 2., "pseudo", grav);

  int i = 0;
  for (; i < nx1; ++i) {
    if (pthermo->GetTemp() < Tmin) break;
    pthermo->GetPrimitive(w1.at(i));

    // set all clouds to zero
    for (int n = 1 + NVAPOR; n < IVX; ++n) w1(n, i) = 0.;

    // move to the next cell
    pthermo->Extrapolate_inplace(dz, "pseudo", grav);
  }

  // Replace adiabatic atmosphere with isothermal atmosphere if temperature
  // is too low
  for (; i < nx1; ++i) {
    pthermo->GetPrimitive(w1.at(i));
    // set all clouds to zero
    for (int n = 1 + NVAPOR; n < IVX; ++n) w1(n, i) = 0.;
    pthermo->Extrapolate_inplace(dz, "isothermal", grav);
  }

  // populate to 3D mesh
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n) {
          phydro->w(n, k, j, i) =
              interp1(pcoord->x1v(i), w1.data() + n * nx1, z1.data(), nx1);
        }
      }

  // add noise
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->w(IVX, k, j, i) = 0.1 * (1. * rand() / RAND_MAX - 0.5);
      }
    }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
