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

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>

Real H2Oratio, CO2ratio, grav;
int iH2O, iH2Oc, iCO2, iCO2c;

Real Ptriple1, Ttriple1;
Real Rd, gammad;
Real x1min, x1max, x2min, x2max;

Real massflux_H2ratio, massflux_CO2ratio;
Real Tm, Ts;

bool fclose(Real x, Real x0) { return std::abs(x - x0) < 1.e-6; }

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1);
  SetUserOutputVariableName(0, "temp");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
      }
}


void BottomInjection(MeshBlock *pmb, Real const time, Real const dt,
                     AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                     AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                     AthenaArray<Real> &s) {
  int is = pmb->is;
  int ie = pmb->ie;

  auto pthermo = Thermodynamics::GetInstance();

  Real p, drhoH2O, drhoH2, drhoCO2;

  Real x1s = pmb->pcoord->x1f(is);

  if (x1s < x1min + pmb->pcoord->dx1f(is)) {
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        // std::cout << pmb->pcoord->x2v(j) << std::endl;
        //  inject at the center of the bottom boundary

        p = pmb->phydro->w(IPR, k, j, is);

        // add water vapor
        drhoH2O =
            dt * 1e-1 /
            pmb->pcoord->dx1f(is);
        // u(iH2O, k, j, is) += drhoH2O;
        // u(IEN, k, j, is) += drhoH2O * (Rd / (gammad - 1.)) *
        //                     pthermo->GetCvRatio(iH2O) * Ttriple1;
        // add dry air (H2)
        drhoH2 = drhoH2O * massflux_H2ratio;
        u(IDN, k, j, is) += drhoH2;
        u(IEN, k, j, is) += drhoH2 * (Rd / (gammad - 1.)) * Ttriple1;

        /* add CO2
        drhoCO2 = drhoH2O * massflux_CO2ratio;
        u(iCO2, k, j, is) += drhoCO2;
        u(IEN, k, j, is) += drhoCO2 * (Rd / (gammad - 1.)) *
                            pthermo->GetCvRatio(iCO2) * Ttriple1;*/
      }
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  BottomInjection(pmb, time, dt, w, r, bcc, u, s);
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  H2Oratio = pin->GetReal("initialcondition", "H2Oratio");
  CO2ratio = pin->GetReal("initialcondition", "CO2ratio");

  grav = -pin->GetReal("hydro", "grav_acc1");

  // index
  iH2O = pthermo->SpeciesIndex("H2O");
  iH2Oc = pthermo->SpeciesIndex("H2O(s)");
  // iCO2 = pthermo->SpeciesIndex("CO2");
  // iCO2c = pthermo->SpeciesIndex("CO2(s)");

  Ptriple1 = pin->GetReal("problem", "Ptriple1");
  Ttriple1 = pin->GetReal("problem", "Ttriple1");
  Rd = pin->GetReal("problem", "Rd");
  gammad = pin->GetReal("hydro", "gamma");
  x1min = pin->GetReal("mesh", "x1min");
  x1max = pin->GetReal("mesh", "x1max");
  x2min = pin->GetReal("mesh", "x2min");
  x2max = pin->GetReal("mesh", "x2max");

  massflux_H2ratio = pin->GetReal("problem", "massflux_H2ratio");
  massflux_CO2ratio = pin->GetReal("problem", "massflux_CO2ratio");

  Tm = pin->GetReal("problem", "Tm");
  Ts = pin->GetReal("problem", "Ts");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();

  // construct 1d atmosphere from bottom up
  std::vector<Real> yfrac(IVX, 0.);
  yfrac[iH2O] = H2Oratio;
  // yfrac[iCO2] = CO2ratio;
  yfrac[0] = 1. - H2Oratio;

  int nx1 = pmy_mesh->mesh_size.nx1;
  Real dz = (x1max - x1min) / (nx1 - 1);
  std::cout << "nx1 = " << nx1 << std::endl;

  AthenaArray<Real> w1, z1;
  w1.NewAthenaArray(NHYDRO, nx1);

  z1.NewAthenaArray(nx1);
  z1(0) = x1min + dz / 2.;
  for (int i = 1; i < nx1; ++i) z1(i) = z1(i - 1) + dz;

  pthermo->SetMassFractions<Real>(yfrac.data());
  pthermo->EquilibrateTP(100., 1.);

  // half a grid to cell center
  pthermo->Extrapolate_inplace(dz / 2., "isothermal", grav);

  for (int i = 0; i < nx1; ++i) {
    pthermo->GetPrimitive(w1.at(i));

    // set all clouds to zero
    for (int n = 1 + NVAPOR; n < IVX; ++n) w1(n, i) = 0.;

    // move to the next cell
    pthermo->Extrapolate_inplace(dz, "isothermal", grav);
  }

  // populate to 3D mesh
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n) {
          phydro->w(n, k, j, i) =
              interp1(pcoord->x1v(i), w1.data() + n * nx1, z1.data(), nx1);
        }
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

}
