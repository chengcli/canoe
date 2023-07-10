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

int iH2O;
Real p0, grav;

struct RootData {
  Real temp_v;
  Variable *qfrac;
};

Real root_func(Real temp, void *aux) {
  auto pthermo = Thermodynamics::GetInstance();

  RootData *pd = static_cast<RootData*>(aux);
  Real& temp_v = pd->temp_v;
  Variable* qfrac = pd->qfrac;

  qfrac->w[IDN] = temp;
  auto rates = pthermo->TryEquilibriumTP(qfrac, iH2O);

  // vapor condensation rate
  qfrac->w[iH2O] += rates[0];

  // cloud concentration rates
  for (int j = 1; j < rates.size(); ++j) {
    qfrac->c[cloud_index_set_[iH2O][j - 1]] += rates[j];

  return thermo->VirtualTemp(qfrac) - temp_v;
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
  auto pthermo = Thermodynamics::GetInstance();

  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      user_out_var(0,j,i) = pthermo->GetTemp(this, ks, j, i);
      user_out_var(1,j,i) = pthermo->PotentialTemp(this, p0, ks, j,i);
      // theta_v
      user_out_var(2,j,i) = user_out_var(1,j,i)*pthermo->RovRd(this, ks, j, i);
      // msv
      user_out_var(3,j,i) = pthermo->MoistStaticEnergy(this, grav*pcoord->x1v(i), ks, j, i);
      // theta_e
      user_out_var(4,j,i) = pthermo->EquivalentPotentialTemp(this, p0, ks, j, i);
      for (int n = 1; n <= NVAPOR; ++n)
        user_out_var(4+n,k,j,i) = pthermo->RelativeHumidity(this, n, ks, j, i);
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  auto pthermo = Thermodynamics::GetInstance();

  Real Ps = p0;
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");
  Real dz = pcoord->dx1f(is);

  // index
  iH2O = pindex->GetVariableId("H2O");
  iH2Oc = pindex->GetVariableId("H2O(c)");

  // index
  auto pindex = IndexMap::GetInstance();
  iH2O = pindex->GetVariableId("H2O");

  Variable qfrac(Variable::Type::MoleFrac);

  // construct a reversible adiabat
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      qfrac.w[iH2O] = qt;
      qfrac.w[IPR] = Ps;
      qfrac.w[IDN] = Ts;
      for (int i = is; i <= ie; ++i) {
        pimpl->DistributeToConserved(qfrac, k, j, i);
        pthermo->Extrapolate(&qfrac, -dz, Thermodynamics::Method::ReversibleAdiabat);
      }
    }

  // add temperature anomaly
  Real temp, Rd = pthermo->GetRd();
  RootData solver;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real L = sqrt(sqr((x2 - xc)/xr) + sqr((x1 - zc)/zr));

        if (L < 1.) {
          pimpl->GatherFromConserved(&qfrac, k, j, i);
          solver.qfrac = &qfrac;
          solver.temp_v = phydro->w(IPR,j,i)/(phydro->w(IDN,j,i)*Rd)*
            (dT*sqr(cos(M_PI*L/2.))/300. + 1.);
          int err = root(qfrac[IDN], qfrac[IDN] + dT, 1.E-8, &temp, root_func, &solver);
          if (err) {
            throw RuntimeError("pgen", "TVSolver doesn't converge");
          } else {
            root_func(temp, &solver);
          }
          pimpl->DistributeToConserved(qfrac, k, j, i);
        }
      }
}
