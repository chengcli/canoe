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

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>

// special includes
#include <special/amars_enroll_vapor_functions_v1.hpp>

Real grav, P0, T0, Tmin;
int iH2O;
int nVars = 7;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(nVars + NVAPOR);
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

  AllocateRealUserMeshBlockDataField(3);
  ruser_meshblock_data[0].NewAthenaArray(ncells2);
  ruser_meshblock_data[1].NewAthenaArray(ncells2);
  ruser_meshblock_data[2].NewAthenaArray(ncells2);

  for (int j = 0; j <= je; ++j) {
    ruser_meshblock_data[0](j) = 0;
    ruser_meshblock_data[1](j) = 230;
    ruser_meshblock_data[2](j) = 0;
  }
}

void MeshBlock::UserWorkInLoop() {
  AthenaArray<Real> &swin = ruser_meshblock_data[0];
  AthenaArray<Real> &ts = ruser_meshblock_data[1];
  // AthenaArray<Real> &ta = ruser_meshblock_data[2];

  AthenaArray<Real> &accumPrecipH2O = ruser_meshblock_data[2];

  double omega = (2 * 3.14159) / 88560;
  double s0 = 1360 * 0.7 * pow(1 / 1.523, 2);
  double alpha_s = 0.3;
  double alpha_a = 0.5;
  double Epsilon_a = 0.5;
  double Epsilon_s = 1;
  double Rho = 1.22;
  double cpAtm = 1005;
  double Sigma = Constants::stefanBoltzmann;
  double Delta_z = 50;
  double cSurf = 100000;
  double time = this->pmy_mesh->time;
  double dt = this->pmy_mesh->dt;

  double precip = 0;
  // double tot_fluxd = 0;
  // for (int i = 0; i<8; ++i){
  //   tot_fluxd += this->prad->
  // }

  for (int j = js; j <= je; ++j) {
    swin(j) = s0 * (1 + std::sin(omega * time));
    // double dTa = ((Epsilon_a * Sigma * (-2 * pow(ta(j), 4) + pow(ts(j), 4)))
    // /
    //               (cpAtm * Delta_z * Rho)) *
    //              dt;
    //  double dTs = ((swin(i) * (1 - alpha_a) * (1 - alpha_s) +
    //                 Epsilon_s * Sigma * (pow(ta(i), 4) - pow(ts(i), 4))) /
    //                cSurf) *
    //               dt;
    double dTs =
        (swin(j) * (1 - alpha_a) * (1 - alpha_s) - Sigma * pow(ts(j), 4)) *
        (dt / cSurf);
    // ta(j) = ta(j) + dTa;
    ts(j) = ts(j) + dTs;

    double iSkim = is + 10;
    for (int i = is; i <= ie; ++i) {
      precip = 0;
      for (int n = 0; n < NCLOUD; ++n) {
        precip = this->pscalars->r(n, ks, j, i);
        this->pscalars->r(n, ks, j, i) = 0;
        if (n == NCLOUD / 2) accumPrecipH2O(j) += precip;
      }
    }
  }
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
