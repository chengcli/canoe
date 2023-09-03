// C++
#include <random>

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

Real phi0, dphi, interval, tseed;
Real omega, vphi, vrad, polarity, skewness;
Real vis;

std::default_random_engine gen;
std::uniform_real_distribution<double> uniform;

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Cubed Sphere");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN, k, j, i) = phi0;
      }
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pexo3 = pimpl->pexo3;
  // Calculate lat, lon, U, V for output and visualization
  for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real Vy = phydro->w(IVY, k, j, i);
        Real Vz = phydro->w(IVZ, k, j, i);
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);
        user_out_var(0, k, j, i) = lat / PI * 180.0;
        user_out_var(1, k, j, i) = lon / PI * 180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = omega * sin(lat);
        Real U, V, U1, V1;
        pexo3->GetUV(&U, &V, Vy, Vz, k, j, i);
        user_out_var(2, k, j, i) = U;
        user_out_var(3, k, j, i) = V;
        // Calculate the vorticity
        Vy = phydro->w(IVY, k + 1, j, i);
        Vz = phydro->w(IVZ, k + 1, j, i);
        pexo3->GetUV(&U, &V, Vy, Vz, k + 1, j, i);
        Vy = phydro->w(IVY, k - 1, j, i);
        Vz = phydro->w(IVZ, k - 1, j, i);
        pexo3->GetUV(&U1, &V1, Vy, Vz, k - 1, j, i);
        Real dz = pcoord->dx3v(k) + pcoord->dx3v(k - 1);
        Real dUdz = (U - U1) / dz;
        Real dVdz = (V - V1) / dz;

        Vy = phydro->w(IVY, k, j + 1, i);
        Vz = phydro->w(IVZ, k, j + 1, i);
        pexo3->GetUV(&U, &V, Vy, Vz, k, j + 1, i);
        Vy = phydro->w(IVY, k, j - 1, i);
        Vz = phydro->w(IVZ, k, j - 1, i);
        pexo3->GetUV(&U1, &V1, Vy, Vz, k, j - 1, i);
        Real dy = pcoord->dx2v(j) + pcoord->dx2v(j - 1);
        Real dUdy = (U - U1) / dy;
        Real dVdy = (V - V1) / dy;

        // Decompose the gradients
        Real dUdlat, dUdlon, dVdlat, dVdlon;
        pexo3->GetUV(&dUdlat, &dUdlon, dUdy, dUdz, k, j, i);
        pexo3->GetUV(&dVdlat, &dVdlon, dVdy, dVdz, k, j, i);

        user_out_var(4, k, j, i) = dVdlon - dUdlat;
        user_out_var(5, k, j, i) =
            (user_out_var(4, k, j, i) + f) / phydro->w(IDN, k, j, i);
      }
}

void MeshBlock::UserWorkInLoop() {
  Real time = pmy_mesh->time;
  Real dt = pmy_mesh->dt;
  auto pexo3 = pimpl->pexo3;
  // check if need to seed a new vortex
  if (time > tseed) {
    Real x1, x2, range;
    Real min_lat = 0.0;
    Real vpol = uniform(gen) > polarity ? 1 : -1;

    if (vpol < 0)  // Anticyclone
      vpol *= skewness;
    Real vlat = -PI / 2.0 + uniform(gen) * PI;
    Real vlon = uniform(gen) * 2. * M_PI;
    for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
      for (int j = js - NGHOST; j <= je + NGHOST; ++j) {
        Real lat, lon, vel;
        pexo3->GetLatLon(&lat, &lon, k, j, is);
        Real radius = pcoord->x1v(is);
        Real dlat = lat - vlat;
        Real dlon = lon - vlon;
        Real a = sin(dlat / 2.0) * sin(dlat / 2.0) +
                 cos(vlat) * cos(lat) * sin(dlon / 2.0) * sin(dlon / 2.0);
        Real c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
        Real a_lon = sin(dlat / 2.0) * sin(dlat / 2.0);
        Real y = 2.0 *
                 atan2(sqrt(a_lon),
                       sqrt(1.0 - a_lon));  // Distance along constant longitude
        Real a_lat = cos(vlat) * cos(lat) * sin(dlon / 2.0) * sin(dlon / 2.0);
        Real x =
            2.0 * atan2(sqrt(a_lat),
                        sqrt(1.0 - a_lat));  // Distance along constant latitude
        // Recover the signs
        if (dlat < 0) y = -y;
        if (dlon < 0) x = -x;

        Real dist = radius * c;
        Real phi = vpol * vphi * exp(-0.5 * (dist * dist) / (vrad * vrad));
        Real fcor = 2. * omega * sin(lat);

        if (vpol > 0)  // cyclone
          vel = phi * dist / (fcor * vrad * vrad);
        else  // anticyclone
          vel = -phi * dist / (fcor * vrad * vrad);
        phydro->w(IDN, k, j, is) += phi;
        // Geostrophic balance
        Real accU = -vel * y / sqrt(x * x + y * y);
        Real accV = vel * x / sqrt(x * x + y * y);
        Real acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, accU, accV, k, j, is);
      }
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord,
                               is - NGHOST, ie + NGHOST, js - NGHOST,
                               je + NGHOST, ks, ke);

    tseed = time - interval * log(uniform(gen));
  }
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

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // forcing parameters
  phi0 = pin->GetReal("problem", "phi0");
  omega = pin->GetReal("problem", "omega");
  interval = pin->GetReal("problem", "interval");
  vrad = pin->GetReal("problem", "vrad");
  vphi = pin->GetReal("problem", "vphi");
  vis = pin->GetOrAddReal("problem", "vis", 0.);
  polarity = pin->GetOrAddReal("problem", "polarity", 0.);
  skewness = pin->GetOrAddReal("problem", "skewness", 1.);

  // next seeding time
  tseed = time;

  // forcing function
  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  // Vort and PV: the derivative is in an unstructured grid, implemented but
  // correctness not verified
  SetUserOutputVariableName(4, "Vort");
  SetUserOutputVariableName(5, "PV");
}
