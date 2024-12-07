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
#include <athena/stride_iterator.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <configure.hpp>
#include <impl.hpp>

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

Real grav, P0, T0, Tmin, prad, hrate, xH2O;
Real Tbot, xtra, ptra, Ttop_tau, Tbot_tau, Trabot_tau;
std::string met;
int iH2O;
bool use_mini, use_mini_ic, use_tra_ic, fix_bot_tra;

//'OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN', 'He'
std::vector<double> vmass = {17.01, 2.02,  18.02, 1.01,  28.01, 44.01, 16.,
                             16.05, 26.04, 17.04, 28.02, 27.03, 4.};
MiniChem *mc;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  if (NVAPOR > 0) {
    AllocateUserOutputVariables(5 + NVAPOR);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "pres");
    SetUserOutputVariableName(3, "thetav");
    SetUserOutputVariableName(4, "mse");
    for (int n = 1; n <= NVAPOR; ++n) {
      std::string name = "rh" + std::to_string(n);
      SetUserOutputVariableName(4 + n, name.c_str());
    }
  } else {
    AllocateUserOutputVariables(3);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "pres");
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, P0, k, j, i);
        user_out_var(2, k, j, i) = phydro->w(IPR, k, j, i);

        if (NVAPOR > 0) {
          // theta_v
          user_out_var(3, k, j, i) =
              user_out_var(1, k, j, i) * pthermo->RovRd(this, k, j, i);
          // mse
          user_out_var(4, k, j, i) =
              pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
          // relative humidity
          for (int n = 1; n <= NVAPOR; ++n)
            user_out_var(4 + n, k, j, i) =
                pthermo->RelativeHumidity(this, n, k, j, i);
        }
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

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);

        Real tem = pthermo->GetTemp(w.at(k,j,i));
        air.ToMoleFraction();
        Real cv =
            get_cv_mass(air, 0, pthermo->GetRd(), pthermo->GetCvRatioMass());

       if (w(IPR, k, j, i) < prad) {
// body cooling 
          du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
                              (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));
// newtonian at the top
           if (tem < Tmin) {
             du(IEN, k, j, i) += dt * (Tmin - tem)/Ttop_tau * w(IDN, k, j, i) * cv;
           }
        }
// relax T at the bottom
         if(i == is) du(IEN, k, j, i) += dt * (Tbot - tem)/Tbot_tau * w(IDN, k, j, i) * cv;
     }
  
  if (NCHEMISTRY > 0) {
    if (use_mini) {
    std::vector<double> vmr(NCHEMISTRY);
    std::vector<double> vmr_sp(NCHEMISTRY - 1);
    std::vector<double> mmr(NCHEMISTRY);
    Real sumVMR;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          // minichem
          Real temp = pthermo->GetTemp(w.at(k, j, i));
          Real pres = w(IPR, k, j, i);

          // change mmr to vmr and normalize
          for (int n = 0; n < NCHEMISTRY; ++n)
            vmr[n] = prim_scalar(NCLOUD + n, k, j, i) / vmass[n];
          sumVMR =
              std::accumulate(vmr.begin(), vmr.end(), static_cast<Real>(0));
          for (auto &value : vmr) value /= sumVMR;

          // call minichem for sp species
          for (int n = 0; n < NCHEMISTRY - 1; ++n) vmr_sp[n] = vmr[n];
          if (temp > 200.) mc->Run(temp, pres, dt, vmr_sp.data(), "NCHO");
          for (int n = 0; n < NCHEMISTRY - 1; ++n) vmr[n] = vmr_sp[n];

          // change vmr to mmr and normalize
          for (int n = 0; n < NCHEMISTRY; ++n) mmr[n] = vmr[n] * vmass[n];
          sumVMR =
              std::accumulate(mmr.begin(), mmr.end(), static_cast<Real>(0));
          for (auto &value : mmr) value /= sumVMR;

          for (int n = 0; n < NCHEMISTRY; ++n) {
            cons_scalar(NCLOUD + n, k, j, i) +=
                phydro->w(IDN, k, j, i) *
                (mmr[n] - prim_scalar(NCLOUD + n, k, j, i));
          }
        }
  } else if (fix_bot_tra) {
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++) {
          if (phydro->w(IPR, k, j, i) > ptra) {
            for (int n = 0; n < NCHEMISTRY; ++n) {
            cons_scalar(NCLOUD + n, k, j, i) +=
                dt/Trabot_tau * phydro->w(IDN, k, j, i) *
                (xtra - prim_scalar(NCLOUD + n, k, j, i));
            }
          }
        }
      }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Tmin = pin->GetReal("problem", "Tmin");
  hrate = pin->GetReal("problem", "hrate") / 86400.;
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate") / 86400.;
  Ttop_tau = pin->GetReal("problem", "Ttop_tau");
  Tbot_tau = pin->GetReal("problem", "Tbot_tau");
  use_mini = pin->GetOrAddBoolean("problem", "use_mini", true);
  use_mini_ic = pin->GetOrAddBoolean("problem", "use_mini_ic", true);
  use_tra_ic = pin->GetOrAddBoolean("problem", "use_tra_ic", true);
  fix_bot_tra = pin->GetOrAddBoolean("problem", "fix_bot_tra", true);
  Trabot_tau = pin->GetReal("problem", "Trabot_tau");
  met = pin->GetOrAddString("problem", "metallicity", "1x");
  xtra = pin->GetReal("problem", "tracer.ppb") / 1.E9;
  ptra = pin->GetReal("problem", "tracer.pre");

  if (NVAPOR > 0) {
    // index
    auto pindex = IndexMap::GetInstance();
    iH2O = pindex->GetVaporId("H2O");
  }
  EnrollUserExplicitSourceFunction(Forcing);

  if (use_mini) {
    // minichem
    mc = new MiniChem();
    mc->SetDataFile("chem_data/mini_chem_data_NCHO.txt");
    mc->SetSpeciesFile("chem_data/mini_chem_sp_NCHO.txt");
    mc->SetNetworkDir("chem_data/" + met + "/");
    mc->SetMetallicityStr(met);
    mc->Initialize();
  }
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
  if (NVAPOR > 0) xH2O = pin->GetReal("problem", "qH2O.ppmv") / 1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    if (NVAPOR > 0) air.w[iH2O] = xH2O;
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
      if (NVAPOR > 0) air.w[iH2O] = xH2O;
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
        //pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav, 1.e-5);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav, 0.);
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

       Tbot = pthermo->GetTemp(this, ks, js, is);;

  if (NCHEMISTRY > 0) {

   if (use_mini_ic) {
    // minichem initialize conserved variables
    // interpolate from ce table
    std::vector<double> vmr_ic(NCHEMISTRY);
    std::vector<double> mmr(NCHEMISTRY);
    std::string ic_file = "chem_data/IC/mini_chem_IC_FastChem_" + met + ".txt";

    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++) {
          double mu;
          double T_in = pthermo->GetTemp(this, k, j, i);
          double P_in = phydro->w(IPR, k, j, i);
          interp_ce_table(NCHEMISTRY, T_in, P_in, vmr_ic.data(), &mu, ic_file);

          // change vmr to mmr and normalize
          for (int n = 0; n < NCHEMISTRY; ++n) mmr[n] = vmr_ic[n] * vmass[n];
          Real sumVMR =
              std::accumulate(mmr.begin(), mmr.end(), static_cast<Real>(0));
          for (auto &value : mmr) value /= sumVMR;

          for (int n = 0; n < NCHEMISTRY; ++n) {
            pscalars->s(NCLOUD + n, k, j, i) +=
                phydro->w(IDN, k, j, i) * mmr[n];
          }
        }
  } else if (use_tra_ic) {
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++) {
          if (phydro->w(IPR, k, j, i) > ptra) {
            for (int n = 0; n < NCHEMISTRY; ++n) {
              pscalars->s(NCLOUD + n, k, j, i) +=
                  phydro->w(IDN, k, j, i) * xtra;
            }
          }
        }
     }
  }
}
