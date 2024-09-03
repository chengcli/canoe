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
// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

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
#include <athena/scalars/scalars.hpp>

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
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>

// minichem
#include <minichem/mini_chem.hpp>

Real grav, P0, T0, Tmin, prad, hrate;
int iH2O;
//'OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN'
std::vector<double> vmass = {17.01, 2.02,  18.02, 1.01,  28.01, 44.01,
                             16.,   16.05, 26.04, 17.04, 28.02, 27.03};
Real mmass = 2.238;  // mean molecular mass in amu
MiniChem *mc;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(3);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "pres");
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
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto pthermo = Thermodynamics::GetInstance();
  auto phydro = pmb->phydro;
  std::vector<Real> vmr(NCHEMISTRY);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);

        air.ToMoleFraction();
        Real cv =
            get_cv_mass(air, 0, pthermo->GetRd(), pthermo->GetCvRatioMass());

        if (w(IPR, k, j, i) < prad) {
          du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
                              (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));

          // minichem
          Real temp = pthermo->GetTemp(w.at(k, j, i));
          Real pres = w(IPR, k, j, i);

          // change mmr to vmr
          for (int n = 0; n < NCHEMISTRY; ++n)
            vmr[n] = prim_scalar(NCLOUD + n, k, j, i) / vmass[n] * mmass;
          //       vmr_from_prim_scalar(vmr, prim_scalar, w, k, j, i);
          // call minichem
          mc->Run(temp, pres, dt, vmr.data(), "NCHO");

          // normalize scale VMR to 1 and change to density
          Real sumVMR =
              std::accumulate(vmr.begin(), vmr.end(), static_cast<Real>(0));
          for (int n = 0; n < NCHEMISTRY; ++n) {
            cons_scalar(NCLOUD + n, k, j, i) =
                phydro->w(IDN, k, j, i) * vmr[n] / sumVMR * vmass[n] / mmass;
          }
        }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");

  Tmin = pin->GetReal("problem", "Tmin");
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate") / 86400.;

  // index
  EnrollUserExplicitSourceFunction(Forcing);

  // minichem
  mc = new MiniChem();
  mc->SetDataFile("chem_data/mini_chem_data_NCHO.txt");
  mc->SetSpeciesFile("chem_data/mini_chem_sp_NCHO.txt");
  mc->SetNetworkDir("chem_data/1x/");
  mc->SetMetallicityStr("1x");
  mc->Initialize();
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

  while (iter++ < max_iter) {
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
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2., "reversible", grav);

      int i = is;
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        if (air.w[IDN] < Tmin) break;
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav, 1.e-5);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        air.w[IVX] = 0.001 * (1. * rand() / RAND_MAX - 0.5);
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "isothermal", grav);
      }
    }

  // minichem initialize conserved variables
  int n_sp = 13;
  // interpolate from ce table
  std::vector<double> vmr_ic(n_sp);
  std::string ic_file = "chem_data/IC/mini_chem_IC_FastChem_1x.txt";

  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
      for (int i = is; i <= ie; i++) {
        double mu;
        double T_in = pthermo->GetTemp(this, k, j, i);
        double P_in = phydro->w(IPR, k, j, i);
        //         std::cout<<"interp  "<<P_in<<" "<<T_in<<std::endl;
        interp_ce_table(n_sp, T_in, P_in, vmr_ic.data(), &mu, ic_file);

        // normalize scale VMR to 1
        Real sumVMR =
            std::accumulate(vmr_ic.begin(), vmr_ic.end(), static_cast<Real>(0));
        for (auto &value : vmr_ic) {
          value /= sumVMR;
        }

        if (NCHEMISTRY > 0) {
          for (int n = 0; n < NCHEMISTRY; ++n) {
            //           std::cout<<"mini  "<<rho<<" "<<n<<" "<<vmr_ic[n]<<"
            //           "<<vmass[n]<<" "<<mu<<std::endl;
            pscalars->s(NCLOUD + n, k, j, i) =
                phydro->w(IDN, k, j, i) * vmr_ic[n] * vmass[n] / mmass;
          }
        }
      }
}
