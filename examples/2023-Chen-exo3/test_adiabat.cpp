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
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>
#include <interface/eos.hpp>

Real Ps, Ts, grav;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto &w = phydro->w;

  auto temp = get_temp(pimpl->peos, w);
  auto pres = get_pres(w);
  auto chi = (get_gammad() - 1.) / get_gammad();
  auto theta = temp * pow(Ps / pres, chi);

  auto temp_a = temp.accessor<Real, 3>();
  auto theta_a = theta.accessor<Real, 3>();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = temp_a[k][j][i];
        user_out_var(1, k, j, i) = theta_a[k][j][i];
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  Ps = pin->GetReal("problem", "Ps");
  Ts = pin->GetReal("problem", "Ts");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  kintera::ThermoX thermo(pimpl->peos->pthermo->options);

  auto temp = Ts * torch::ones({ncells3, ncells2}, torch::kDouble);
  auto pres = Ps * torch::ones({ncells3, ncells2}, torch::kDouble);
  auto xfrac = torch::ones({ncells3, ncells2, 1}, torch::kDouble);

  // half a grid to cell center
  thermo->extrapolate_ad(temp, pres, xfrac, grav, pcoord->dx1f(is) / 2.);

  auto w = get_all(phydro->w);
  for (int i = is; i <= ie; ++i) {
    auto conc = thermo->compute("TPX->V", {temp, pres, xfrac});
    w[IPR].select(2, i) = pres;
    w[IDN].select(2, i) = thermo->compute("V->D", {conc});
    w.slice(0, 1, IVX).select(3, i) = thermo->compute("X->Y", {xfrac});
    thermo->extrapolate_ad(temp, pres, xfrac, grav, pcoord->dx1f(i));
  }

  // Change primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
