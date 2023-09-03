// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Cubed Sphere");

  auto pexo3 = pimpl->pexo3;

  // Input field is a uniform 10 m/s zonal wind
  Real V = 0.0;
  Real U = 0.0;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real R = pcoord->x1v(i);
        Real R0 = 5E5;
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        Real rad = (PI / 2.0 - lat) * R;
        if ((rad < R0) && (lat > PI / 4.0)) {
          phydro->w(IDN, k, j, i) = 500.0 + 10.0 * cos(PI / 2.0 * rad / R0);
        } else {
          phydro->w(IDN, k, j, i) = 500.0;
        }
        rad = (PI / 2.0 + lat) * R;
        if ((rad < R0) && (lat < -PI / 4.0)) {
          phydro->w(IDN, k, j, i) = 500.0 + 10.0 * cos(PI / 2.0 * rad / R0);
        }

        Real Vy, Vz;
        pexo3->GetVyVz(&Vy, &Vz, U, V, k, j, i);
        phydro->w(IVY, k, j, i) = Vy;
        phydro->w(IVZ, k, j, i) = Vz;
      }
  app->Log("Done with problem generator");
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  int is = pmb->is;
  Real omega = 2. * PI / (24. * 3600.);
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real cF2, cF3;
        pexo3->CalculateCoriolisForce2(j, k, w(IVY, k, j, i), w(IVZ, k, j, i),
                                       omega, w(IDN, k, j, i), &cF2, &cF3);
        u(IVY, k, j, i) += dt * cF2;
        u(IVZ, k, j, i) += dt * cF3;
      }
    }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real C = sqrt(1.0 + x * x);
        Real D = sqrt(1.0 + y * y);
        Real delta = 1.0 / (1.0 + x * x + y * y);

        Real Vy = phydro->w(IVY, k, j, i);
        Real Vz = phydro->w(IVZ, k, j, i);
        Real lat, lon;
        pimpl->pexo3->GetLatLon(&lat, &lon, k, j, i);
        user_out_var(0, k, j, i) = lat / PI * 180.0;
        user_out_var(1, k, j, i) = lon / PI * 180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 0.0;
        Real U, V;
        pimpl->pexo3->GetUV(&U, &V, Vy, Vz, k, j, i);
        user_out_var(2, k, j, i) = U;
        user_out_var(3, k, j, i) = V;
        user_out_var(4, k, j, i) = delta / (C * D);
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  SetUserOutputVariableName(4, "sqrtg");
}

// Uncomment the following to turn on Coriolis force
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Real gamma = pin->GetReal("hydro", "gamma");
  EnrollUserExplicitSourceFunction(Forcing);
}
