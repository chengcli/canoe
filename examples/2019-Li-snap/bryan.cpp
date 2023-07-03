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
#include <athena/parameter_input.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <impl.hpp>

// climath
#include <climath/core.h>
#include <climath/root.h>
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include "../utils/utils.hpp"

int iH2O, iH2Oc;
Real p0, grav;

struct RootData {
  Real temp_v;
  Real beta;
  Real delta;
  Real p3;
  Real t3;
  Real eps;
  Real *qfrac;
};

Real root_func(Real temp, void *aux) {
  RootData *pd = static_cast<RootData*>(aux);
  Real& temp_v = pd->temp_v;
  Real& beta = pd->beta;
  Real& delta = pd->delta;
  Real& p3 = pd->p3;
  Real& t3 = pd->t3;
  Real& eps = pd->eps;
  Real *qfrac = pd->qfrac;

  qfrac[IDN] = temp;
  Real rate = GasCloudIdeal(qfrac, iH2O, iH2Oc, t3, p3, 0., beta, delta);
  qfrac[iH2O] -= rate;
  qfrac[iH2Oc] += rate;
  Real f1 = 1. - qfrac[iH2Oc];
  Real f2 = 1. + (qfrac[iH2O] + qfrac[iH2Oc])*(eps - 1.);
  return temp*f1/f2 - temp_v;
};

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "theta_v");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "theta_e");
  SetUserOutputVariableName(5, "rh");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  auto pthermo = pimpl->pthermo;

  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      user_out_var(0,j,i) = pthermo->GetTemp(phydro->w.at(j,i));
      user_out_var(1,j,i) = pthermo->PotentialTemp(phydro->w.at(j,i), p0);
      // theta_v
      user_out_var(2,j,i) = user_out_var(1,j,i)*pthermo->RovRd(phydro->w.at(j,i));
      // msv
      user_out_var(3,j,i) = pthermo->MoistStaticEnergy(phydro->w.at(j,i), grav*pcoord->x1v(i));
      // theta_e
      user_out_var(4,j,i) = pthermo->EquivalentPotentialTemp(phydro->w.at(j,i), p0);
      for (int n = 1; n <= NVAPOR; ++n)
        user_out_var(4+n,k,j,i) = pthermo->RelativeHumidity(phydro->w.at(k,j,i), n);
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");

  // index
  iH2O = pindex->GetVariableId("H2O");
  iH2Oc = pindex->GetVariableId("H2O(c)");

  RootData solver;
  solver.p3 = pin->GetReal("thermodynamics", "Ptriple1");
  solver.t3 = pin->GetReal("thermodynamics", "Ttriple1");
  solver.eps = pthermo->GetMassRatio(iH2Oc);
  solver.beta = pthermo->GetBeta(iH2Oc);
  solver.delta = pthermo->GetDelta(iH2Oc);

  Real Rd = pthermo->GetRd();
  Real gamma = peos->GetGamma();
  Real cpd = gamma/(gamma - 1.)*Rd;
  
  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  Real dz = (x1max - x1min)/(nx1 - 1);
  Real **w1, *z1;
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];

  w1[0][iH2O] = qt;
  Real Ps = p0;
  pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::reversible);
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  auto pcloud = pimpl->cloudq[0];

  for (int i = is; i <= ie; ++i) {
    Real buf[Variable::Size];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, Variable::Size);
    for (int j = js; j <= je; ++j) {
      for (int n = 0; n < NHYDRO; ++n)
        phydro->w(n,j,i) = buf[n];
      pcloud->r(iH2Oc - NHYDRO,j,i) = buf[iH2Oc];
      phydro->w(IVX,j,i) = 0.;
      phydro->w(IVY,j,i) = 0.;
    }
  }

  // add temperature anomaly
  Real qfrac[Variable::Size]; Real temp;
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(sqr((x2 - xc)/xr) + sqr((x1 - zc)/zr));
      if (L < 1.) {
        pimpl->ToMolarFraction(qfrac, ks, j, i);
        solver.qfrac = qfrac;
        solver.temp_v = phydro->w(IPR,j,i)/(phydro->w(IDN,j,i)*Rd)*
          (dT*sqr(cos(M_PI*L/2.))/300. + 1.);
        int err = root(qfrac[IDN], qfrac[IDN] + dT, 1.E-8, &temp, root_func, &solver);
        if (err) {
          throw RuntimeError("pgen", "TVSolver doesn't converge");
        } else {
          root_func(temp, &solver);
        }
        pimpl->FromMolarFraction(qfrac, ks, j, i);
      }
    }
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
}
