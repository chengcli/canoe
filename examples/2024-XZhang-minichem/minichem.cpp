// C++ headers
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

// yaml
#include <yaml-cpp/yaml.h>

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

// snap
#include <snap/eos/ideal_moist.hpp>
#include <snap/mesh/meshblock.hpp>

// canoe
#include <configure.h>

#include <impl.hpp>
#include <interface/eos.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// minichem
#include <minichem/mini_chem.hpp>

Real grav, P0, T0, Tmin, prad, hrate, xH2O;
Real Tbot, xtra, ptra, Ttop_tau, Tbot_tau, Trabot_tau;
std::string met;
bool use_mini, use_mini_ic, use_tra_ic, fix_bot_tra;

//'OH','H2','H2O','H','CO','CO2','O','CH4','C2H2','NH3','N2','HCN', 'He'
std::vector<double> vmass = {17.01, 2.02,  18.02, 1.01,  28.01, 44.01, 16.,
                             16.05, 26.04, 17.04, 28.02, 27.03, 4.};
MiniChem *mc;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  int nvapor = 0;

  if (nvapor > 0) {
    AllocateUserOutputVariables(5 + nvapor);
    SetUserOutputVariableName(0, "temp");
    SetUserOutputVariableName(1, "theta");
    SetUserOutputVariableName(2, "pres");
    SetUserOutputVariableName(3, "thetav");
    SetUserOutputVariableName(4, "mse");
    for (int n = 1; n <= nvapor; ++n) {
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
  auto w = phydro->w;
  auto temp = get_temp(pimpl->peos, w);
  auto pres = get_pres(w);
  auto chi = (get_gammad() - 1.) / get_gammad();
  auto theta = temp * pow(P0 / pres, chi);

  auto temp_a = temp.accessor<Real, 3>();
  auto theta_a = theta.accessor<Real, 3>();
  auto pres_a = pres.accessor<Real, 3>();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = temp_a[k][j][i];
        user_out_var(1, k, j, i) = theta_a[k][j][i];
        user_out_var(2, k, j, i) = pres_a[k][j][i];

        /*if (NVAPOR > 0) {
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
        }*/
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &du,
             AthenaArray<Real> &cons_scalar) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  auto phydro = pmb->phydro;

  auto temp = get_temp(pmb->pimpl->peos, w);
  auto cv = get_cv(pmb->pimpl->peos, w);

  auto temp_a = temp.accessor<Real, 3>();
  auto cv_a = cv.accessor<Real, 3>();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        auto tem = temp_a[k][j][i];
        auto cv = cv_a[k][j][i];
        if (w(IPR, k, j, i) < prad) {
          // body cooling
          du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
                              (1. + 1.E-4 * sin(2. * M_PI * rand() / RAND_MAX));
          // newtonian at the top
          if (tem < Tmin) {
            du(IEN, k, j, i) +=
                dt * (Tmin - tem) / Ttop_tau * w(IDN, k, j, i) * cv;
          }
        }
        // relax T at the bottom
        if (i == is)
          du(IEN, k, j, i) +=
              dt * (Tbot - tem) / Tbot_tau * w(IDN, k, j, i) * cv;
      }

  if (NCHEM > 0) {
    if (use_mini) {
      std::vector<double> vmr(NCHEM);
      std::vector<double> vmr_sp(NCHEM - 1);
      std::vector<double> mmr(NCHEM);
      Real sumVMR;
      auto temp_a = temp.accessor<Real, 3>();

      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) {
            // minichem
            Real temp = temp_a[k][j][i];
            Real pres = w(IPR, k, j, i);

            // change mmr to vmr and normalize
            for (int n = 0; n < NCHEM; ++n)
              vmr[n] = prim_scalar(n, k, j, i) / vmass[n];
            sumVMR =
                std::accumulate(vmr.begin(), vmr.end(), static_cast<Real>(0));
            for (auto &value : vmr) value /= sumVMR;

            // call minichem for sp species
            for (int n = 0; n < NCHEM - 1; ++n) vmr_sp[n] = vmr[n];
            if (temp > 200.) mc->Run(temp, pres, dt, vmr_sp.data(), "NCHO");
            for (int n = 0; n < NCHEM - 1; ++n) vmr[n] = vmr_sp[n];

            // change vmr to mmr and normalize
            for (int n = 0; n < NCHEM; ++n) mmr[n] = vmr[n] * vmass[n];
            sumVMR =
                std::accumulate(mmr.begin(), mmr.end(), static_cast<Real>(0));
            for (auto &value : mmr) value /= sumVMR;

            for (int n = 0; n < NCHEM; ++n) {
              cons_scalar(n, k, j, i) +=
                  phydro->w(IDN, k, j, i) * (mmr[n] - prim_scalar(n, k, j, i));
            }
          }
    } else if (fix_bot_tra) {
      for (int k = ks; k <= ke; k++)
        for (int j = js; j <= je; j++)
          for (int i = is; i <= ie; i++) {
            if (phydro->w(IPR, k, j, i) > ptra) {
              for (int n = 0; n < NCHEM; ++n) {
                cons_scalar(n, k, j, i) += dt / Trabot_tau *
                                           phydro->w(IDN, k, j, i) *
                                           (xtra - prim_scalar(n, k, j, i));
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
  Tbot = pin->GetOrAddReal("problem", "Tbot", 0.);
  use_mini = pin->GetOrAddBoolean("problem", "use_mini", true);
  use_mini_ic = pin->GetOrAddBoolean("problem", "use_mini_ic", true);
  use_tra_ic = pin->GetOrAddBoolean("problem", "use_tra_ic", true);
  fix_bot_tra = pin->GetOrAddBoolean("problem", "fix_bot_tra", true);
  Trabot_tau = pin->GetReal("problem", "Trabot_tau");
  met = pin->GetOrAddString("problem", "metallicity", "1x");
  xtra = pin->GetReal("problem", "tracer.ppb") / 1.E9;
  ptra = pin->GetReal("problem", "tracer.pre");

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

torch::Tensor setup_moist_adiabatic_profile(std::string infile) {
  auto config = YAML::LoadFile(infile);
  auto Ps = config["problem"]["Ps"].as<double>(1.e5);
  auto Ts = config["problem"]["Ts"].as<double>(300.);
  auto Tmin = config["problem"]["Tmin"].as<double>(200.);
  auto grav = -config["forcing"]["const-gravity"]["grav1"].as<double>();

  // initialize the block
  auto block = snap::MeshBlock(snap::MeshBlockOptions::from_yaml(infile));

  // useful modules
  auto phydro = block->phydro;
  auto pcoord = phydro->pcoord;
  auto peos = phydro->peos;
  auto m = block->named_modules()["hydro.eos.thermo"];
  auto thermo_y = std::dynamic_pointer_cast<kintera::ThermoYImpl>(m);

  // dimensions and indices
  int nc3 = pcoord->x3v.size(0);
  int nc2 = pcoord->x2v.size(0);
  int nc1 = pcoord->x1v.size(0);
  int ny = thermo_y->options.species().size() - 1;
  int nvar = peos->nvar();

  // construct an adiabatic atmosphere
  kintera::ThermoX thermo_x(thermo_y->options);

  auto temp = Ts * torch::ones({nc3, nc2}, torch::kFloat64);
  auto pres = Ps * torch::ones({nc3, nc2}, torch::kFloat64);
  auto xfrac = torch::zeros({nc3, nc2, 1 + ny}, torch::kFloat64);
  auto w = torch::zeros({nvar, nc3, nc2, nc1}, torch::kFloat64);

  // read in compositions
  for (int i = 1; i <= ny; ++i) {
    auto name = thermo_y->options.species()[i];
    auto xmixr = config["problem"]["x" + name].as<double>(0.);
    xfrac.select(2, i) = xmixr;
  }

  // dry air mole fraction
  xfrac.select(2, 0) = 1. - xfrac.narrow(-1, 1, ny).sum(-1);

  // adiabatic extrapolate half a grid to cell center
  int is = pcoord->is();
  int ie = pcoord->ie();
  auto dz = pcoord->dx1f[is].item<double>();
  thermo_x->extrapolate_ad(temp, pres, xfrac, grav, dz / 2.);

  int i = is;
  int nvapor = thermo_x->options.vapor_ids().size();
  int ncloud = thermo_x->options.cloud_ids().size();
  for (; i <= ie; ++i) {
    auto conc = thermo_x->compute("TPX->V", {temp, pres, xfrac});

    w[IPR].select(2, i) = pres;
    w[IDN].select(2, i) = thermo_x->compute("V->D", {conc});

    auto result = thermo_x->compute("X->Y", {xfrac});
    w.narrow(0, snap::ICY, ny).select(3, i) =
        thermo_x->compute("X->Y", {xfrac});

    if ((temp < Tmin).any().item<double>()) break;
    dz = pcoord->dx1f[i].item<double>();
    thermo_x->extrapolate_ad(temp, pres, xfrac, grav, dz);
  }

  // isothermal extrapolation
  for (; i <= ie; ++i) {
    auto mu = (thermo_x->mu * xfrac).sum(-1);
    dz = pcoord->dx1f[i].item<double>();
    pres *= exp(-grav * mu * dz / (kintera::constants::Rgas * temp));
    auto conc = thermo_x->compute("TPX->V", {temp, pres, xfrac});
    w[IPR].select(2, i) = pres;
    w[IDN].select(2, i) = thermo_x->compute("V->D", {conc});
    w.narrow(0, snap::ICY, ny).select(3, i) =
        thermo_x->compute("X->Y", {xfrac});
  }

  // add noise
  w[IVX] += 0.01 * torch::rand_like(w[IVX]);
  w[IVY] += 0.01 * torch::rand_like(w[IVY]);

  return w;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // auto infile = pin->GetString("problem", "config_file");
  // auto w = setup_moist_adiabatic_profile(infile)

  srand(Globals::my_rank + time(0));

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pin->GetReal("hydro", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;
  std::cout << "Rd = " << Rd << ", cp = " << cp << std::endl;

  // set up an adiabatic atmosphere
  int max_iter = 400, iter = 0;
  Real Ttol = pin->GetOrAddReal("problem", "init_Ttol", 0.01);

  // estimate surface temperature and pressure
  auto x1min = pcoord->x1f(is);
  auto x1max = pcoord->x1f(ie + 1);

  Real Ts = T0 - grav / cp * x1min;
  Real Ps = P0 * pow(Ts / T0, cp / Rd);

  // Loop over the grids and set initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real temp = T0 - grav * x1 / cp;
        phydro->w(IPR, k, j, i) = P0 * pow(temp / T0, cp / Rd);
        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * temp);
        phydro->w(IVX, k, j, i) = 0.001 * (1. * rand() / RAND_MAX - 0.5);
      }

  auto temp = get_temp(pimpl->peos, phydro->w);
  auto temp_a = temp.accessor<Real, 3>();

  if (std::abs(Tbot - temp_a[ks][js][is]) > 1.e-2) {
    std::cout << "Please set Tbot = " << temp_a[ks][js][is]
              << " K in the input file" << std::endl;
    throw std::runtime_error("Wrong Tbot");
  }

  std::cout << "ptra = " << ptra << " Pa" << std::endl;

  if (NCHEM > 0) {
    if (use_mini_ic) {
      // minichem initialize conserved variables
      // interpolate from ce table
      std::vector<double> vmr_ic(NCHEM);
      std::vector<double> mmr(NCHEM);
      std::string ic_file =
          "chem_data/IC/mini_chem_IC_FastChem_" + met + ".txt";

      for (int k = ks; k <= ke; k++)
        for (int j = js; j <= je; j++)
          for (int i = is; i <= ie; i++) {
            double mu;
            double T_in = temp_a[k][j][i];
            double P_in = phydro->w(IPR, k, j, i);
            interp_ce_table(NCHEM, T_in, P_in, vmr_ic.data(), &mu, ic_file);

            // change vmr to mmr and normalize
            for (int n = 0; n < NCHEM; ++n) mmr[n] = vmr_ic[n] * vmass[n];
            Real sumVMR =
                std::accumulate(mmr.begin(), mmr.end(), static_cast<Real>(0));
            for (auto &value : mmr) value /= sumVMR;

            for (int n = 0; n < NCHEM; ++n) {
              pscalars->s(n, k, j, i) += phydro->w(IDN, k, j, i) * mmr[n];
            }
          }
    } else if (use_tra_ic) {
      for (int k = ks; k <= ke; k++)
        for (int j = js; j <= je; j++)
          for (int i = is; i <= ie; i++) {
            if (phydro->w(IPR, k, j, i) > ptra) {
              for (int n = 0; n < NCHEM; ++n) {
                pscalars->s(n, k, j, i) += phydro->w(IDN, k, j, i) * xtra;
              }
            }
          }
    }
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
