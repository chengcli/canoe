// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/core.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

//! variables used in initial condition
Real grav, P0, T0, gamma;
Real xCO2, xCH4, xN2O;
Real zt, zq1, zq2;
Real qt, q0;
int iH2O = 1, iCO2 = 0, iCH4 = 1, iN2O = 2, iO3 = 3;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "rh_H2O");

  AllocateRealUserMeshBlockDataField(1);

  // CO2, CH4, N2O, O3 in mole fraction
  ruser_meshblock_data[0].NewAthenaArray(4, ncells1);
  for (int i = is; i <= ie; ++i) {
    ruser_meshblock_data[0](iCO2, i) = xCO2;
    ruser_meshblock_data[0](iCH4, i) = xCH4;
    ruser_meshblock_data[0](iN2O, i) = xN2O;
  }
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
        user_out_var(4, k, j, i) = pthermo->RelativeHumidity(this, 1, k, j, i);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
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
  gamma = pin->GetReal("problem", "gamma");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  Real Rd = pthermo->GetRd();
  Real Tv0 = T0 * (1. + 0.608 * q0);
  Real Tvt = Tv0 - gamma * zt;
  Real Pt = P0 * pow(Tvt / Tv0, grav / (Rd * gamma));

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
              P0 * pow((Tv0 - gamma * z) / Tv0, grav / (Rd * gamma));
          Tv = Tv0 - gamma * z;
        } else {
          phydro->w(iH2O, k, j, i) = qt;
          phydro->w(IPR, k, j, i) = Pt * exp(-grav * (z - zt) / (Rd * Tvt));
          Tv = Tvt;
        }

        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * Tv);

        // pa -> hpa
        Real pre = phydro->w(IPR, k, j, i) / 100.;
        // ppmv -> mole fraction
        ruser_meshblock_data[0](iO3, i) =
            1.e-6 * g1 * pow(pre, g2) * exp(-pre / g3);
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
