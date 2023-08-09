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
Real flux_h2o, flux_dry, flux_heat, flux_width;

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
        AirParcel air(AirParcel::Type::MassFrac);
        pimpl->GatherFromPrimitive(&air, k, j, i);
        user_out_var(6, k, j, i) = air.w[iH2O] + air.c[iH2Oc];
      }
}

double sat_vapor_p_H2O(double T) {
  double betal = 24.845, gammal = -2.1735, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

// water svp
void Thermodynamics::enrollVaporFunctionH2O() {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1_[iH2O][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_H2O(qfrac.w[IDN]);
  };
}

void CondensateGravity(MeshBlock *pmb, Real const time, Real const dt,
                       AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                       AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                       AthenaArray<Real> &s) {
  // acceleration in 1-direction
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
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

  // surface plume source
  Real x2mid =
      (pmb->pmy_mesh->mesh_size.x2min + pmb->pmy_mesh->mesh_size.x2max) / 2.;
  Real x2min = x2mid - flux_width / 2.;
  Real x2max = x2mid + flux_width / 2.;

  Real x3mid =
      (pmb->pmy_mesh->mesh_size.x3min + pmb->pmy_mesh->mesh_size.x3max) / 2.;
  Real x3min = x3mid - flux_width / 2.;
  Real x3max = x3mid + flux_width / 2.;

  int is = pmb->is;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real x3 = pmb->pcoord->x3v(k);
      Real x2 = pmb->pcoord->x2v(j);

      bool source_x2 =
          (!pmb->pmy_mesh->f2) || (x2 > x2min) && (x2 < x2max) ? true : false;
      bool source_x3 =
          (!pmb->pmy_mesh->f3) || (x3 > x3min) && (x3 < x3max) ? true : false;

      if (source_x2 && source_x3) {
        u(iH2O, k, j, is) += dt * flux_h2o / pmb->pcoord->dx1f(is);
        u(IDN, k, j, is) += dt * flux_dry / pmb->pcoord->dx1f(is);
        u(IEN, k, j, is) += dt * flux_heat / pmb->pcoord->dx1f(is);
      }
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");

  flux_h2o = pin->GetReal("problem", "flux_h2o");
  flux_dry = pin->GetReal("problem", "flux_dry");
  flux_heat = pin->GetReal("problem", "flux_heat");
  flux_width = pin->GetReal("problem", "flux_width");

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

  AirParcel air(AirParcel::Type::MassFrac);
  air.w[iH2O] = qt;
  air.c[iH2Oc] = 0.;

  air.ToMoleFraction();
  qt = air.w[iH2O];

  // construct a reversible adiabat
  srand(Globals::my_rank + time(0));
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.w[iH2O] = qt;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;
      air.c[iH2Oc] = 0.;

      // half a grid to cell center
      pthermo->Extrapolate(&air, pcoord->dx1f(is) / 2.,
                           Thermodynamics::Method::PseudoAdiabat, grav);

      for (int i = is; i <= ie; ++i) {
        // add noise
        air.w[IVX] = 0.01 * (1. * rand() / RAND_MAX - 0.5);
        pimpl->DistributeToConserved(air, k, j, i);
        pthermo->Extrapolate(&air, pcoord->dx1f(i),
                             Thermodynamics::Method::PseudoAdiabat, grav);
      }
    }
}
