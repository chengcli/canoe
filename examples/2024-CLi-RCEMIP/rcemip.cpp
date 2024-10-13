// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// climath
#include <climath/core.h>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

//! variables used in initial condition
Real grav, P0, T0, gamma_tv;
Real xCO2, xCH4, xN2O;
Real zt, zq1, zq2;
Real qt, q0;
int iH2O, iH2Oc, iH2Op;
int iCO2 = 0, iCH4 = 1, iN2O = 2, iO3 = 3;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(11);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "tempv");
  SetUserOutputVariableName(2, "enthalpy");
  SetUserOutputVariableName(3, "entropy");
  SetUserOutputVariableName(4, "intEng");

  SetUserOutputVariableName(5, "theta");
  SetUserOutputVariableName(6, "thetav");
  SetUserOutputVariableName(7, "mse");

  SetUserOutputVariableName(8, "rh_H2O");
  SetUserOutputVariableName(9, "theta_e");
  SetUserOutputVariableName(10, "qtol");

  // CO2, CH4, N2O, O3 in mole fraction
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(4, ncells1);
  for (int i = is; i <= ie; ++i) {
    ruser_meshblock_data[0](iCO2, i) = xCO2;
    ruser_meshblock_data[0](iCH4, i) = xCH4;
    ruser_meshblock_data[0](iN2O, i) = xN2O;
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

        user_out_var(5, k, j, i) = potential_temp(pthermo, w.at(k, j, i), P0);
        // theta_v
        user_out_var(6, k, j, i) =
            user_out_var(5, k, j, i) * pthermo->RovRd(w.at(k, j, i));

        // mse
        user_out_var(7, k, j, i) =
            moist_static_energy(pthermo, w.at(k, j, i), grav * pcoord->x1v(i));

        // rh
        user_out_var(8, k, j, i) = relative_humidity(pthermo, w.at(k, j, i))[1];
        // theta_e
        user_out_var(9, k, j, i) = equivalent_potential_temp(
            pthermo, w.at(k, j, i), user_out_var(8, k, j, i), P0);

        // total mixing ratio
        user_out_var(10, k, j, i) =
            w(iH2O, k, j, i) + w(iH2Oc, k, j, i) + w(iH2Op, k, j, i);
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  // heat capacity
  Real cv = 717.;           // J/kg/K
  Real dTdt = 2. / 86400.;  // K/s
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        u(IEN, k, j, i) -= cv * w(IDN, k, j, i) * dTdt * dt;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  grav = -pin->GetReal("hydro", "grav_acc1");

  xCO2 = pin->GetReal("problem", "xCO2");
  xCH4 = pin->GetReal("problem", "xCH4");
  xN2O = pin->GetReal("problem", "xN2O");

  zt = pin->GetReal("problem", "zt");
  zq1 = pin->GetReal("problem", "zq1");
  zq2 = pin->GetReal("problem", "zq2");

  qt = pin->GetReal("problem", "qt");
  q0 = pin->GetReal("problem", "q0");

  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  gamma_tv = pin->GetReal("problem", "gamma_tv");

  // index
  iH2O = pthermo->SpeciesIndex("H2O");
  iH2Oc = pthermo->SpeciesIndex("H2O(l)");
  iH2Op = pthermo->SpeciesIndex("H2O(p)");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  srand(Globals::my_rank + time(0));

  Application::Logger app("main");
  app->Log("ProblemGenerator: rcemip");

  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  Real Rd = pthermo->GetRd();
  Real Tv0 = T0 * (1. + 0.608 * q0);
  Real Tvt = Tv0 - gamma_tv * zt;
  Real Pt = P0 * pow(Tvt / Tv0, grav / (Rd * gamma_tv));

  // O3 parameters
  Real g1 = 3.6478;  // ppmv hPa−g2
  Real g2 = 0.83209;
  Real g3 = 11.3515;  // hPa

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real z = pcoord->x1v(i);
        Real Tv;

        if (z < zt) {
          phydro->w(iH2O, k, j, i) = q0 * exp(-z / zq1) * exp(-sqr(z / zq2));
          phydro->w(IPR, k, j, i) =
              P0 * pow((Tv0 - gamma_tv * z) / Tv0, grav / (Rd * gamma_tv));
          Tv = Tv0 - gamma_tv * z;
        } else {
          phydro->w(iH2O, k, j, i) = qt;
          phydro->w(IPR, k, j, i) = Pt * exp(-grav * (z - zt) / (Rd * Tvt));
          Tv = Tvt;
        }

        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * Tv);

        // random noise
        phydro->w(IVX, k, j, i) = 0.1 * (1. * rand() / RAND_MAX - 0.5);

        // pa -> hpa
        Real pre = phydro->w(IPR, k, j, i) / 100.;

        // ppmv -> mole fraction
        ruser_meshblock_data[0](iO3, i) =
            1.e-6 * g1 * pow(pre, g2) * exp(-pre / g3);
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
