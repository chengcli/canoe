#include <fstream>

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

// microphysics
#include <microphysics/microphysics.hpp>

// snap
// #include <snap/thermodynamics/calc_surf_evaporation_rates.cpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// special includes
#include <special/amars_enroll_vapor_functions_v1.hpp>

Real grav, P0, T0, Tmin;
int iH2O;
int nNewVars = 7;
int nTotVars = nNewVars + 4;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(nTotVars + NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3 + n, name.c_str());
  }
  SetUserOutputVariableName(4 + NVAPOR, "btemp");
  SetUserOutputVariableName(4 + NVAPOR + 1, "accumPrecipH2O");
  // for now, amd is a scalar and associated with one phase. later, make amd[n]
  // for each phase
  SetUserOutputVariableName(4 + NVAPOR + 2, "lH2Oamd");
  SetUserOutputVariableName(4 + NVAPOR + 3, "sH2Oamd");
  SetUserOutputVariableName(4 + NVAPOR + 4, "lH2Ogel");
  SetUserOutputVariableName(4 + NVAPOR + 5, "sH2Ogel");

  AllocateRealUserMeshBlockDataField(nNewVars);
  ruser_meshblock_data[0].NewAthenaArray(ncells2);
  ruser_meshblock_data[1].NewAthenaArray(ncells2);
  ruser_meshblock_data[2].NewAthenaArray(ncells2);
  ruser_meshblock_data[3].NewAthenaArray(ncells2);
  ruser_meshblock_data[4].NewAthenaArray(ncells2);
  ruser_meshblock_data[5].NewAthenaArray(ncells2);
  ruser_meshblock_data[6].NewAthenaArray(ncells2);

  for (int j = 0; j <= je; ++j) {
    ruser_meshblock_data[0](j) = 0;
    ruser_meshblock_data[1](j) = 230;
    ruser_meshblock_data[2](j) = 0;
    ruser_meshblock_data[3](j) = 0;
    ruser_meshblock_data[4](j) = 100;
    ruser_meshblock_data[5](j) = 0;
    ruser_meshblock_data[6](j) = 0;
  }
}

void MeshBlock::UserWorkInLoop() {
  AthenaArray<Real> &swin = ruser_meshblock_data[0];
  AthenaArray<Real> &ts = ruser_meshblock_data[1];

  AthenaArray<Real> &accumPrecipH2O = ruser_meshblock_data[2];
  AthenaArray<Real> &lH2Oamd = ruser_meshblock_data[3];
  AthenaArray<Real> &sH2Oamd = ruser_meshblock_data[4];
  AthenaArray<Real> &lH2Ogel = ruser_meshblock_data[5];
  AthenaArray<Real> &sH2Ogel = ruser_meshblock_data[6];

  double omega = (2 * 3.14159) / 88560;
  double s0 = 1360 * 0.7 * pow(1 / 1.523, 2);
  double alpha_s = 0.3;
  double alpha_a = 0.5;
  double Sigma = Constants::stefanBoltzmann;
  double cSurf = 100000;
  double time = this->pmy_mesh->time;
  double dt = this->pmy_mesh->dt;

  double dTs = 0;
  double precip = 0;
  double tot_fluxdn = 0;
  int numBands = this->pimpl->prad->GetNumBands();

  bool H2OisLiquid;
  auto pthermo = Thermodynamics::GetInstance();
  double rholH2O = 1000;
  double rhosH2O = 910;
  double dz = pcoord->x1f(is + 1) - pcoord->x1f(is);
  int iSkim;

  for (int j = js; j <= je; ++j) {
    swin(j) = s0 * (1 + std::sin(omega * time));

    // get the flux from each band and add it up
    tot_fluxdn = 0;
    for (int n = 0; n < numBands; ++n) {
      tot_fluxdn += this->pimpl->prad->GetBand(n)->bflxdn(ks, j, is);
    }

    dTs = (swin(j) * (1 - alpha_a) * (1 - alpha_s) + tot_fluxdn -
           Sigma * pow(ts(j), 4)) *
          (dt / cSurf);
    ts(j) = ts(j) + dTs;

    // track precip
    if (ts(j) > 273)
      H2OisLiquid = true;
    else
      H2OisLiquid = false;

    // only skim off the first layer of precip
    iSkim = is;
    for (int i = is; i <= iSkim; ++i) {
      precip = 0;
      for (int n = NCLOUD / 2; n < NCLOUD; ++n) {
        precip = this->pscalars->r(n, ks, j, i);
        this->pscalars->r(n, ks, j, i) = 0;
        if (n == 1) {
          if (H2OisLiquid)
            lH2Oamd(j) += precip * this->phydro->w(IDN, ks, j, is) * dz;
          else
            sH2Oamd(j) += precip * this->phydro->w(IDN, ks, j, is) * dz;
          accumPrecipH2O(j) += precip;
        }
      }
    }

    // get the current air parcel
    AirParcel air(AirParcel::Type::MoleFrac);
#pragma omp simd
    for (int n = 0; n < NHYDRO; ++n) air.w[n] = this->phydro->w(n, ks, j, is);
#pragma omp simd
    for (int n = 0; n < NCLOUD; ++n)
      air.c[n] = this->pimpl->pmicro->u(n, ks, j, is);

    Real Mbar = pthermo->GetMu(this, ks, j, is);
    // make sure air has T in IDN slot
    air.w[IDN] = (air.w[IPR] * Mbar) / (air.w[IDN] * Constants::Rgas);

    for (int i = 1; i <= NVAPOR; ++i) {
      std::vector<Real> rates(1 + pthermo->GetCloudIndexSet(i).size(), 0.);
      //(cmetz) check this value of CDE, see Hartmann pg 117
      rates = pthermo->CalcSurfEvapRates(air, i, lH2Oamd(j), ts(j), dTs, cSurf,
                                         dt, 3e-3, Mbar);
      std::cout << "liquid rates: " << rates[0] << std::endl;
      lH2Oamd(j) += -rates[0];
      this->phydro->w(i, ks, j, is) += rates[0] * (Mbar / pthermo->GetMu(i)) /
                                       (this->phydro->w(IDN, ks, j, is) * dz);
      if (this->phydro->w(i, ks, j, is) < 0) this->phydro->w(i, ks, j, is) = 0;

      rates = pthermo->CalcSurfEvapRates(air, i, sH2Oamd(j), ts(j), dTs, cSurf,
                                         dt, 3e-3, Mbar);
      std::cout << "solid rates: " << rates[0] << std::endl;
      sH2Oamd(js) += -rates[0];
      this->phydro->w(i, ks, j, is) += rates[0] * (Mbar / pthermo->GetMu(i)) /
                                       (this->phydro->w(IDN, ks, j, is) * dz);
      if (this->phydro->w(i, ks, j, is) < 0) this->phydro->w(i, ks, j, is) = 0;
    }

    lH2Ogel(j) = lH2Oamd(j) / rholH2O;
    sH2Ogel(j) = sH2Oamd(j) / rhosH2O;

  }  // end j loop
}

// void Mesh::UserWorkAfterLoop(ParameterInput *pin) {

// for (int b=0; b<nblocal; ++b) {
//     MeshBlock *pmb = my_blocks(b);

//  AthenaArray<Real> &accumPrecipH2O = pmb->ruser_meshblock_data[2];
//  AthenaArray<Real> &accumPrecipH2S = pmb->ruser_meshblock_data[3];
//  AthenaArray<Real> &accumPrecipSO2 = pmb->ruser_meshblock_data[4];

//    BoundaryValues *pbval = pmb->pbval;
//    int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je,
//        kl = pmb->ks, ku = pmb->ke;

//  double precip = 0;
//    for (int n = 0; n < NCLOUD; ++n){
//      for (int j = 0; j <= ju; ++j) {
//        precip = pmb->pscalars->r(n, kl, j, il);
//        pmb->pscalars->r(n, kl, j, il) = 0;

//        if(n=3) accumPrecipH2O(j) += precip;
//        if(n=4) accumPrecipH2S(j) += precip;
//        if(n=5) accumPrecipSO2(j) += precip;
//      }
//  }

//}

//}

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
        user_out_var(4 + NVAPOR, k, j, i) = ruser_meshblock_data[1](j);
        user_out_var(4 + NVAPOR + 1, k, j, i) = ruser_meshblock_data[2](j);

        user_out_var(4 + NVAPOR + 2, k, j, i) = ruser_meshblock_data[3](j);
        user_out_var(4 + NVAPOR + 3, k, j, i) = ruser_meshblock_data[4](j);
        user_out_var(4 + NVAPOR + 4, k, j, i) = ruser_meshblock_data[5](j);
        user_out_var(4 + NVAPOR + 5, k, j, i) = ruser_meshblock_data[6](j);
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
        auto &&air = AirParcelHelper::gather_from_primitive(pmb, k, j, i);

        // Real cv = pthermo->GetCvMass(air, 0);

        // if (w(IPR, k, j, i) < prad) {
        //   du(IEN, k, j, i) += dt * hrate * w(IDN, k, j, i) * cv *
        //                       (1. + 1.E-4 * sin(2. * M_PI * rand() /
        //                       RAND_MAX));
        // }

        // if (air.w[IDN] < Tmin) {
        //   u(IEN,k,j,i) += w(IDN,k,j,i)*cv*(Tmin - temp)/sponge_tau*dt;
        // }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");

  Tmin = pin->GetReal("problem", "Tmin");

  // index
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVaporId("H2O");
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: amarsw_crm");

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
  // Real xCO2 = pin->GetReal("problem", "qCO2.ppmv") / 1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    air.w[iH2O] = xH2O;
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
      // air.w[iCO2] = xCO2;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2., "pseudo", grav);

      int i = is;
      for (; i <= ie; ++i) {
        if (air.w[IDN] < Tmin) break;
        air.w[IVX] = 0.1 * sin(2. * M_PI * rand() / RAND_MAX);
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav, 1.e-3);
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_conserved(this, k, j, i, air);
        pthermo->Extrapolate(&air, pcoord->dx1f(i), "isothermal", grav);
      }

      peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b, phydro->w,
                                 pfield->bcc, pcoord, is, ie, j, j, k, k);

      pimpl->prad->CalFlux(this, k, j, is, ie + 1);
    }
}
