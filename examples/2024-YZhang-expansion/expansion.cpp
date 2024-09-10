// @sect3{Include files}

// These input files are just like those in the @ref straka, so additional
// comments are not required.
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <impl.hpp>

// snap
#include <snap/stride_iterator.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

// @sect3{Preamble}
#define sqr(x) ((x) * (x))

// We only need one global variables here, the surface pressure
Real p0;

// Same as that in @ref straka, make outputs of temperature and potential
// temperature.
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

// Set temperature and potential temperature.
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  auto pthermo = Thermodynamics::GetInstance();
  auto &w = phydro->w;

  Real gamma = peos->GetGamma();
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(this, p0, k, j, i);
      }
}

// Initialize surface pressure from input file.
void Mesh::InitUserMeshData(ParameterInput *pin) {
  p0 = pin->GetReal("problem", "p0");
}

// @sect3{Initial condition}

// We do not need forcings other than gravity in this problem,
// so we go directly to the initial condition.
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Similar to @ref straka, read variables in the input file
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Ts = pin->GetReal("problem", "Ts");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;

  Real xc = pin->GetReal("problem", "xc");
  Real yc = pin->GetReal("problem", "yc");
  Real zc = pin->GetReal("problem", "zc");
  Real a = pin->GetReal("problem", "a");
  Real dT = pin->GetReal("problem", "dT");
  Real dP = pin->GetReal("problem", "dP");

  // Loop over the grids and set initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real temp = Ts;
        phydro->w(IPR, k, j, i) = p0;
        for (int n = 0; n < 5; ++n) {
          Real x = 0.04 * rand() / RAND_MAX - 0.02;
          Real y = 0.04 * rand() / RAND_MAX - 0.02;
          Real z = 0.04 * rand() / RAND_MAX - 0.02;
          Real r = sqrt(sqr(x3 - y - yc) + sqr(x2 - x - xc) + sqr(x1 - z - zc));
          if (r <= a) {
            temp = dT;
            phydro->w(IPR, k, j, i) = dP;
          }
        }
        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * temp);
        phydro->w(IVX, k, j, i) = phydro->w(IVY, k, j, i) = 0.;
      }

  // Change primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
