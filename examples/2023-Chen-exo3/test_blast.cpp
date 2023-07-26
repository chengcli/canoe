
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>
#include <configure.hpp>
#include <exo3/cubed_sphere.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Input field is a uniform 10 m/s zonal wind
  Real V = 0.0;
  Real U = 0.0;
  CubedSphere cs(this);

  std::cout << "======Start of Problem generator======" << loc.lx2 << "||"
            << loc.lx3 << std::endl;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x = tan(pcoord->x2v(j));
        Real y = tan(pcoord->x3v(k));
        Real R = pcoord->x1v(i);
        Real R0 = 5E5;
        Real lat, lon;
        cs.GetLatLon(&lat, &lon, k, j, i);
        Real rad = (PI / 2.0 - lat) * R;
        // if ((rad < R0) && (lat>PI/4.0))
        //   phydro->w(IDN,k,j,i) = 500.0 + 10.0 * cos(PI/2.0*rad/R0);
        // else
        phydro->w(IDN, k, j, i) = 1.0;  // / (R*R);
        Real Vy, Vz;
        cs.GetVyVz(&Vy, &Vz, U, V, k, j, i);
        phydro->w(IVY, k, j, i) = Vy;
        phydro->w(IVZ, k, j, i) = Vz;
      }
  std::cout << "Done with problem generator..." << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Calculate lat, lon, U, V for output and visualization
  CubedSphere cs(this);
  for (int k = ks - NGHOST; k <= ke + NGHOST; ++k)
    for (int j = js - NGHOST; j <= je + NGHOST; ++j)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        Real Vy = phydro->w(IVY, k, j, i);
        Real Vz = phydro->w(IVZ, k, j, i);
        Real lat, lon;
        cs.GetLatLon(&lat, &lon, k, j, i);
        user_out_var(0, k, j, i) = lat / PI * 180.0;
        user_out_var(1, k, j, i) = lon / PI * 180.0;
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 0.0;
        Real U, V;
        cs.GetUV(&U, &V, Vy, Vz, k, j, i);
        user_out_var(2, k, j, i) = U;
        user_out_var(3, k, j, i) = V;
        user_out_var(4, k, j, i) = pcoord->GetFace1Area(k, j, i);
      }
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "lat");
  SetUserOutputVariableName(1, "lon");
  SetUserOutputVariableName(2, "U");
  SetUserOutputVariableName(3, "V");
  SetUserOutputVariableName(4, "Area");
}
