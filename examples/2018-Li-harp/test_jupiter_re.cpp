// C/C++
#include <ctime>
#include <sstream>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>
#include <athena/stride_iterator.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// math
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// global parameters
Real grav, P0, T0, Z0, Tmin;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real dq[1 + NVAPOR], rh;
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) =
            pthermo->GetTemp(this, k, j, i);
        user_out_var(1, k, j, i) =
            pthermo->PotentialTemp(this, P0, k, j, i);
      }
}

// w, primitive hydro
// r, primitive scalar
// u, conserved hydro
// s, conserved scalar
void Forcing(MeshBlock *pmb, Real const time, Real const dt,
             AthenaArray<Real> const &w, AthenaArray<Real> const &r,
             AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
             AthenaArray<Real> &s) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // damp kinetic energy
        u(IM1, k, j, i) -= w(IDN, k, j, i) * w(IM1, k, j, i) / 5.;
        u(IM2, k, j, i) -= w(IDN, k, j, i) * w(IM2, k, j, i) / 5.;
        u(IM3, k, j, i) -= w(IDN, k, j, i) * w(IM3, k, j, i) / 5.;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  grav = -pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Z0 = pin->GetOrAddReal("problem", "Z0", 0.);
  Tmin = pin->GetReal("hydro", "min_tem");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Application::Logger app("main");

  app->Log("ProblemGenerator: test_jupiter_re");
  Real gamma = pin->GetReal("hydro", "gamma");

  // construct a 1D adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  Real **w1, *z1, *p1, *t1;
  Real dz = (x1max - x1min) / pmy_mesh->mesh_size.nx1 / 2.;
  size_t nx1 = 2 * pmy_mesh->mesh_size.nx1 + 1;
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real[nx1];
  p1 = new Real[nx1];
  t1 = new Real[nx1];

  // estimate surface temperature and pressure
  auto pthermo = Thermodynamics::GetInstance();

  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;
  Real Ts = T0 - grav / cp * (x1min - Z0);
  Real Ps = P0 * pow(Ts / T0, cp / Rd);

  // set up an adiabatic atmosphere
  int max_iter = 200, iter = 0;
  Real dlnp = pcoord->dx1f(is)/H0;

  Variable var(Variable::Type::MoleFrac);

  // estimate surface temperature and pressure
  Real Ps = P0*exp(-x1min/H0);
  Real Ts = T0*pow(Ps/P0, Rd/cp);
  Real xH2O = pin->GetReal("problem", "qH2O.ppmv")/1.E6;
  Real xNH3 = pin->GetReal("problem", "qNH3.ppmv")/1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    var.w[iH2O] = xH2O;
    var.w[iNH3] = xNH3;
    var.w[IPR] = Ps;
    var.w[IDN] = Ts;

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate(&var, -dlnp/2., Thermodynamics::Method::DryAdiabat);
      if (var.w[IPR] < P0) break;
    }

    // extrapolate down to where var is
    pthermo->Extrapolate(&var, log(P0/var.w[IPR]),
        Thermodynamics::Method::DryAdiabat);

    // make up for the difference
    Ts += T0 - var.w[IDN];
    if (std::abs(T0 - var.w[IDN]) < 0.01) break;

    app->Log("Iteration #", iter);
    app->Log("T", var.w[IDN]);
  }

  if (iter > max_iter) {
    throw RuntimeError("ProblemGenerator", "maximum iteration reached");
  }

  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      var.SetZero();
      var.w[iH2O] = xH2O;
      var.w[iNH3] = xNH3;
      var.w[IPR] = Ps;
      var.w[IDN] = Ts;

      int i = is;
      for (; i <= ie; ++i) {
        pimpl->DistributeToPrimitive(var, k, j, i);
        pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::DryAdiabat);
        if (var.w[IDN] < Tmin) break;
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
      pthermo->Extrapolate(&var, dlnp, Thermodynamics::Method::DryAdiabat);
      for (; i <= ie; ++i) {
        pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::Isothermal);
        pimpl->DistributeToPrimitive(var, k, j, i);
      }
    }

  // if requested
  // modify atmospheric stability
  Real adlnTdlnP = pin->GetOrAddReal("problem", "adlnTdlnP", 0.);

  if (adlnTdlnP != 0.) {
    Real pmin = pin->GetOrAddReal("problem", "adlnTdlnP.pmin", 1.);
    Real pmax = pin->GetOrAddReal("problem", "adlnTdlnP.pmax", 20.);

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        int ibegin = find_pressure_level_lesser(pmax, phydro->w, k, j, is, ie);
        int iend = find_pressure_level_lesser(pmin, phydro->w, k, j, is, ie);

        pimpl->GatherFromPrimitive(&var, k, j, ibegin);

        for (int i = ibegin; i < iend; ++i) {
          pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::DryAdiabat,
              0., adlnTdlnP);
          pimpl->DistributeToPrimitive(var, k, j, i+1);
        }
      }
  }

  // add noise
  unsigned int seed = time(NULL);
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        phydro->w(IVX, k, j, i) = 0.01 * (1. * rand_r(&seed) / RAND_MAX - 0.5);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
