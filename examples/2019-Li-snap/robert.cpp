/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2021
 * Contact: chengcli@umich.edu
 * Reference: Robert et al., 1992
 * -------------------------------------------------------------------------------------
 */

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
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(w.at(k, j, i), p0);
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
  auto pthermo = Thermodynamics::GetInstance();

  // Similar to @ref straka, read variables in the input file
  Real gamma = pin->GetReal("hydro", "gamma");
  Real grav = -phydro->hsrc.GetG1();
  Real Ts = pin->GetReal("problem", "Ts");
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;

  Real xc = pin->GetReal("problem", "xc");
  Real yc = pin->GetReal("problem", "yc");
  Real zc = pin->GetReal("problem", "zc");
  Real s = pin->GetReal("problem", "s");
  Real a = pin->GetReal("problem", "a");
  Real dT = pin->GetReal("problem", "dT");

  // Whether to do a uniform bubble or not.
  bool uniform_bubble =
      pin->GetOrAddBoolean("problem", "uniform_bubble", false);

  // Loop over the grids and set initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real r = sqrt((x3 - yc) * (x3 - yc) + (x2 - xc) * (x2 - xc) +
                      (x1 - zc) * (x1 - zc));
        Real temp = Ts - grav * x1 / cp;
        phydro->w(IPR, k, j, i) = p0 * pow(temp / Ts, cp / Rd);
        if (r <= a)
          temp += dT * pow(phydro->w(IPR, k, j, i) / p0, Rd / cp);
        else if (!uniform_bubble)
          temp += dT * exp(-(r - a) * (r - a) / (s * s)) *
                  pow(phydro->w(IPR, k, j, i) / p0, Rd / cp);
        phydro->w(IDN, k, j, i) = phydro->w(IPR, k, j, i) / (Rd * temp);
        phydro->w(IVX, k, j, i) = phydro->w(IVY, k, j, i) = 0.;
      }

  // Change primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
