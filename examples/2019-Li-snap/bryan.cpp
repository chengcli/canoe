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
#include <snap/thermodynamics/thermodynamics.hpp>

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

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, p0, k, j, i);
        // theta_v
        user_out_var(2, k, j, i) =
            user_out_var(1, j, i) * pthermo->RovRd(this, k, j, i);
        // mse
        user_out_var(3, k, j, i) =
            pthermo->MoistStaticEnergy(this, grav * pcoord->x1v(i), k, j, i);
        // theta_e
        user_out_var(4, k, j, i) =
            pthermo->EquivalentPotentialTemp(this, p0, iH2O, k, j, i);
        user_out_var(5, k, j, i) =
            pthermo->RelativeHumidity(this, iH2O, k, j, i);
        // total mixing ratio
        Variable air(Variable::Type::MassFrac);
        pimpl->GatherFromPrimitive(&air, k, j, i);
        user_out_var(6, k, j, i) = air.w[iH2O] + air.c[iH2Oc];
      }
}

double sat_vapor_p_H2O(double T) {
  double betal = 24.845, gammal = -2.1735, tr = 273.16, pr = 611.7;
  return svph2o(T / tr, pr, betal, gammal);
}

// water svp
void Thermodynamics::enrollVaporFunctionH2O() {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1_[iH2O][0] = [](Variable const &qfrac, int, int) {
    return sat_vapor_p_H2O(qfrac.w[IDN]);
  };
}

void CondensateGravity(MeshBlock *pmb, Real const time, Real const dt,
                       AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                       AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                       AthenaArray<Real> &s) {
  // acceleration in 1-direction
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real qd = 1.;
#pragma omp simd reduction(+ : qd)
        for (int n = 1; n <= NVAPOR; ++n) qd += -w(n, k, j, i);

        Real rho_dry = w(IDN, k, j, i) * qd;
        Real rho_cloud = 0.;
#pragma omp simd
        for (int n = 0; n < NCLOUD; ++n) rho_cloud += rho_dry * r(n, k, j, i);

        Real src = -dt * rho_cloud * grav;
        u(IM1, k, j, i) += src;
        if (NON_BAROTROPIC_EOS) u(IEN, k, j, i) += src * w(IVX, k, j, i);
      }
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");
  EnrollUserExplicitSourceFunction(CondensateGravity);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pindex = IndexMap::GetInstance();

  Real Ps = p0;
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");

  // index
  iH2O = pindex->GetVaporId("H2O");
  iH2Oc = pindex->GetCloudId("H2O(c)");

  Variable air(Variable::Type::MassFrac);
  air.w[iH2O] = qt;
  air.c[iH2Oc] = 0.;

  air.ToMoleFraction();
  qt = air.w[iH2O];

  // construct a reversible adiabat
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[iH2O] = qt;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;
      air.c[iH2Oc] = 0.;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2.,
                           Thermodynamics::Method::ReversibleAdiabat, grav);

      for (int i = is; i <= ie; ++i) {
        pimpl->DistributeToConserved(air, k, j, i);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::ReversibleAdiabat, grav);
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
          pimpl->GatherFromConserved(&air, k, j, i);
          Real temp_v = air.w[IDN] * pthermo->RovRd(air);
          temp_v *= 1. + dT * sqr(cos(M_PI * L / 2.)) / 300.;

          Real temp;
          int err = root(air.w[IDN], air.w[IDN] + dT, 1.E-8, &temp,
                         [&pthermo, &air, temp_v](Real temp) {
                           air.w[IDN] = temp;

                           auto rates =
                               pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);
                           air.w[iH2O] += rates[0];
                           air.c[iH2Oc] += rates[1];

                           return temp * pthermo->RovRd(air) - temp_v;
                         });

          if (err) throw RuntimeError("pgen", "TVSolver doesn't converge");

          air.w[IDN] = temp;
          auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);
          air.w[iH2O] += rates[0];
          air.c[iH2Oc] += rates[1];

          pimpl->DistributeToConserved(air, k, j, i);
        }
      }
}
