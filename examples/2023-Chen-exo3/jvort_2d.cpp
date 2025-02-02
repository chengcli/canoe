// athena
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <random>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exo3
#include <exo3/cubed_sphere.hpp>
#include <exo3/cubed_sphere_utility.hpp>

Real phi0, dphi, Li, omega, change_tau, radius, cool_tau;
int M; // number of forcing components

AthenaArray<Real> As, Bs, Lats, Lons; // Amplitudes, phases, latitudes, longitudes
std::default_random_engine gen;
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
              AthenaArray<Real> const &w, const AthenaArray<Real> &prim_scalar,
              AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
              AthenaArray<Real> &cons_scalar) {
  auto pexo3 = pmb->pimpl->pexo3;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lat, lon;
        pexo3->GetLatLon(&lat, &lon, k, j, i);

        // Coriolis force
        Real f = 2. * omega * sin(lat);
        Real U, V;
        pexo3->GetUV(&U, &V, w(IVY, k, j, i), w(IVZ, k, j, i), k, j, i);
        Real ll_acc_U = f * V;
        Real ll_acc_V = -f * U;
        Real acc1, acc2, acc3;
        pexo3->GetVyVz(&acc2, &acc3, ll_acc_U, ll_acc_V, k, j, i);
        pexo3->ContravariantVectorToCovariant(j, k, acc2, acc3, &acc2, &acc3);
        u(IM2, k, j, i) += dt * w(IDN, k, j, i) * acc2;
        u(IM3, k, j, i) += dt * w(IDN, k, j, i) * acc3;

        // Forcing
        Real forcing = 0.0;
        Real wn = 2.0 * M_PI * radius / Li;
        for (int m = 0; m < M; ++m) {
          Real arg = acos(sin(lat)*sin(Lats(m)) + cos(lat)*cos(Lats(m))*cos(lon - Lons(m)));
          forcing += As(m) * cos(wn * arg + Bs(m)) * dphi;
        }
        u(IDN, k, j, i) += dt * forcing;
      }
      
  // Evolve As, Bs, Lats, Lons
  for (int m = 0; m < M; ++m) {
    As(m) += dt * normal(gen) / change_tau * 1.0;
    Bs(m) += dt * uniform(gen) / change_tau * 2.0 * M_PI;
    Real theta = uniform(gen) * 2.0 * M_PI;  // random direction
    Real dist = fabs(normal(gen)) * Li * dt / change_tau;  // random distance
    
    // Calculate new latitude using spherical trig
    Real phi = Lats(m);
    Real phi_new = asin(sin(phi)*cos(dist/radius) + 
                       cos(phi)*sin(dist/radius)*cos(theta));
    
    // Calculate new longitude using spherical trig
    Real lambda = Lons(m);
    Real denom = cos(dist/radius) - sin(phi)*sin(phi_new);
    Real dlambda = atan2(sin(theta)*sin(dist/radius)*cos(phi), denom);
    Real lambda_new = lambda + dlambda;
    // If latitude exceeds poles, flip to opposite hemisphere and reverse longitude
    if (phi_new > M_PI/2) {
        phi_new = M_PI - phi_new;
        lambda_new += M_PI;
    } else if (phi_new < -M_PI/2) {
        phi_new = -M_PI - phi_new; 
        lambda_new += M_PI;
    }
    
    // Constrain longitude between -pi and pi
    while (lambda_new > M_PI) lambda_new -= 2*M_PI;
    while (lambda_new < -M_PI) lambda_new += 2*M_PI;
    Lats(m) = phi_new;
    Lons(m) = lambda_new;
  }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");
  app->Log("ProblemGenerator: Uniform field, cubed sphere, phi0 = ", phi0);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN, k, j, i) = phi0;
      }
  app->Log("Done with problem generator");
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
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
        Real forcing = 0.0;
        Real wn = 2.0 * M_PI * radius / Li;
        for (int m = 0; m < M; ++m) {
          Real arg = acos(sin(lat)*sin(Lats(m)) + cos(lat)*cos(Lats(m))*cos(lon - Lons(m)));
          forcing += As(m) * cos(wn * arg + Bs(m)) * dphi;
        }
        user_out_var(5, k, j, i) = forcing;
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  SetUserOutputVariableName(4, "sqrtg");
  SetUserOutputVariableName(5, "forcing");
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // forcing parameters
  phi0 = pin->GetReal("problem", "phi0");
  dphi = pin->GetReal("problem", "dphi");
  Li = pin->GetReal("problem", "Li");
  M = pin->GetInteger("problem", "M");
  omega = pin->GetReal("problem", "omega");
  change_tau = pin->GetReal("problem", "change_tau");
  radius = pin->GetReal("problem", "radius");
  cool_tau = pin->GetReal("problem", "cool_tau");

  As.NewAthenaArray(M);
  Bs.NewAthenaArray(M);
  Lats.NewAthenaArray(M);
  Lons.NewAthenaArray(M);
  for (int m = 0; m < M; ++m) {
    As(m) = 1.0;
    Bs(m) = uniform(gen) * 2.0 * M_PI;
    // Use fibonacci sphere algorithm for even distribution
    Real golden_ratio = (1 + sqrt(5)) / 2;
    Lats(m) = acos(1 - 2*(m+0.5)/M);  // theta
    Lons(m) = 2*M_PI * m/golden_ratio; // phi
  }
  EnrollUserExplicitSourceFunction(Forcing);
}
