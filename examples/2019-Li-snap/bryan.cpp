/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2023
 * Contact: chengcli@umich.edu
 * Reference: Bryan and Fritsch, 2002
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

#include <climath/root.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// special includes
// #include "bryan_vapor_functions.hpp"

int iH2O, iH2Oc;
Real p0, grav;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "theta_v");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "theta_e");
  SetUserOutputVariableName(5, "rh");
  SetUserOutputVariableName(6, "qtol");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(w.at(k, j, i), p0);
        // theta_v
        user_out_var(2, k, j, i) =
            user_out_var(1, j, i) * pthermo->RovRd(w.at(k, j, i));
        /* mse
        user_out_var(3, k, j, i) =
            pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
        // theta_e
        user_out_var(4, k, j, i) =
            pthermo->EquivalentPotentialTemp(this, p0, iH2O, k, j, i);
        user_out_var(5, k, j, i) =
            pthermo->RelativeHumidity(this, iH2O, k, j, i);*/
        // total mixing ratio
        user_out_var(6, k, j, i) = w(iH2O, k, j, i) + w(iH2Oc, k, j, i);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pindex = IndexMap::GetInstance();

  grav = -pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");

  // index
  iH2O = pindex->GetVaporId("H2O");
  iH2Oc = pindex->GetCloudId("H2O(c)");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  Real Ps = p0;
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");

  // construct a reversible adiabat
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      w(IDN, k, j, is) = Ts;
      w(IPR, k, j, is) = Ps;
      w(iH2O, k, j, is) = qt;
      w(iH2Oc, k, j, is) = 0.;

      pthermo->SetPrimitive(w.at(k, j, is));
      // half a grid to cell center
      pthermo->Extrapolate_inplace(pcoord->dx1f(is) / 2., "reversible", grav);

      for (int i = is; i <= ie; ++i) {
        pthermo->GetPrimitive(w.at(k, j, i));
        pthermo->Extrapolate_inplace(pcoord->dx1f(i), "reversible", grav);
      }
    }

  // add temperature anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real L = sqrt(sqr((x2 - xc) / xr) + sqr((x1 - zc) / zr));

        if (L < 1.) {
          pthermo->SetPrimitive(w.at(k, j, i));
          Real temp = pthermo->GetTemp();
          Real temp_v = temp * pthermo->RovRd();
          temp_v *= 1. + dT * sqr(cos(M_PI * L / 2.)) / 300.;

          int err = root(temp, temp + dT, 1.E-8, &temp,
                         [&pthermo, temp_v](Real temp) {
                           pthermo->SetTemperature(temp);
                           pthermo->EquilibrateTP();
                           Real rovrd = pthermo->RovRd();
                           return temp * rovrd - temp_v;
                         });

          if (err) throw RuntimeError("pgen", "TVSolver doesn't converge");

          pthermo->SetTemperature(temp);
          pthermo->EquilibrateTP();
          pthermo->GetPrimitive(w.at(k, j, i));
        }
      }

  peos->PrimitiveToConserved(w, pfield->bcc, phydro->u, pcoord, is, ie, js, je,
                             ks, ke);
}
