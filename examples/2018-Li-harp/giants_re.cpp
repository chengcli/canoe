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

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <dirty.hpp>
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
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) = pthermo->GetTemp(w.at(k, j, i));
        user_out_var(1, k, j, i) = pthermo->PotentialTemp(w.at(k, j, i), P0);
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
  app->Log("ProblemGenerator: jupiter_re");

  auto pthermo = Thermodynamics::GetInstance();

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real H0 = pcoord->GetPressureScaleHeight();

  // request temperature and pressure
  app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;

  // index
  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");
  int iNH3 = pindex->GetVaporId("NH3");

  // set up an adiabatic atmosphere
  int max_iter = 200, iter = 0;
  Real dlnp = pcoord->dx1f(is) / H0;

  AirParcel air(AirParcel::Type::MoleFrac);

  // estimate surface temperature and pressure
  Real Ps = P0 * exp(-x1min / H0);
  Real Ts = T0 * pow(Ps / P0, Rd / cp);
  Real xH2O = pin->GetReal("problem", "qH2O.ppmv") / 1.E6;
  Real xNH3 = pin->GetReal("problem", "qNH3.ppmv") / 1.E6;

  while (iter++ < max_iter) {
    // read in vapors
    air.w[iH2O] = xH2O;
    air.w[iNH3] = xNH3;
    air.w[IPR] = Ps;
    air.w[IDN] = Ts;

    // stop at just above P0
    for (int i = is; i <= ie; ++i) {
      pthermo->Extrapolate(&air, -dlnp / 2., "dry");
      if (air.w[IPR] < P0) break;
    }

    // extrapolate down to where var is
    pthermo->Extrapolate(&air, log(P0 / air.w[IPR]), "dry");

    // make up for the difference
    Ts += T0 - air.w[IDN];
    if (std::abs(T0 - air.w[IDN]) < 0.01) break;

    app->Log("Iteration #", iter);
    app->Log("T", air.w[IDN]);
  }

  if (iter > max_iter) {
    throw RuntimeError("ProblemGenerator", "maximum iteration reached");
  }

  // construct atmosphere from bottom up
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      air.SetZero();
      air.w[iH2O] = xH2O;
      air.w[iNH3] = xNH3;
      air.w[IPR] = Ps;
      air.w[IDN] = Ts;

      int i = is;
      for (; i <= ie; ++i) {
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);

        pthermo->Extrapolate(&air, -dlnp, "dry");
        if (air.w[IDN] < Tmin) break;
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      pthermo->Extrapolate(&air, dlnp, "dry");
      for (; i <= ie; ++i) {
        pthermo->Extrapolate(&air, -dlnp, "isothermal");
        AirParcelHelper::distribute_to_primitive(this, k, j, i, air);
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

        auto &&air = AirParcelHelper::gather_from_primitive(this, k, j, ibegin);
        air.ToMoleFraction();

        for (int i = ibegin; i < iend; ++i) {
          pthermo->Extrapolate(&air, -dlnp, "dry", 0., adlnTdlnP);
          AirParcelHelper::distribute_to_primitive(this, k, j, i + 1, air);
        }
      }
  }

  // set chemical tracers
  auto ptracer = pimpl->ptracer;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
