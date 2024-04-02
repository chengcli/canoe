// C/C++
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp> 
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// set up an adiabatic atmosphere
void construct_atmosphere(MeshBlock *pmb, ParameterInput *pin, Real NH3ppmv, Real T0) {
  Application::Logger app("main");
  // app->Log("ProblemGenerator: juno");

  app->Log("NH3.ppmv", NH3ppmv);
  app->Log("T0", T0);

  auto pmy_mesh= pmb->pmy_mesh;
  auto pthermo = Thermodynamics::GetInstance();
  auto pcoord = pmb->pcoord;

  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js, ks = pmb->ks;
  int je = pmb->je, ke = pmb->ke;

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  Real H0 = pcoord->GetPressureScaleHeight();

  Real P0 = pin->GetReal("mesh", "ReferencePressure");

  Real Tmin = pin->GetReal("problem", "Tmin");
  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;

  // index
  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");
  int iNH3 = pindex->GetVaporId("NH3");

  // app->Log("index of H2O", iH2O);
  // app->Log("index of NH3", iNH3);

  // set up an adiabatic atmosphere
  int max_iter = 200, iter = 0;
  Real dlnp = pcoord->dx1f(is) / H0;

  AirParcel air(AirParcel::Type::MoleFrac);

  // estimate surface temperature and pressure
  Real Ps = P0 * exp(-x1min / H0);
  Real Ts = T0 * pow(Ps / P0, Rd / cp);
  Real xH2O = pin->GetReal("problem", "qH2O.ppmv") / 1.E6;
//   Real xNH3 = pin->GetReal("problem", "qNH3.ppmv") / 1.E6;
  Real xNH3 = NH3ppmv/ 1.E6;
  // app->Log("xH2O", xH2O);
  // app->Log("xNH3", xNH3);
  // app->Log("x1min", x1min);
  // app->Log("P0", P0);
  // app->Log("Rd", Rd);
  // app->Log("gamma", gamma);

  Real rh_max_nh3 = pin->GetOrAddReal("problem", "rh_max.NH3", 1.);

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

    // extrapolate down to where air is
    pthermo->Extrapolate(&air, log(P0 / air.w[IPR]), "dry");

    // make up for the difference
    Ts += T0 - air.w[IDN];
    if (std::abs(T0 - air.w[IDN]) < 0.01) break;

    // app->Log("Iteration #", iter);
    // app->Log("T", air.w[IDN]);
  }

  if (iter > max_iter) {
    throw RuntimeError("ProblemGenerator", "maximum iteration reached");
  }

  // construct atmosphere from bottom up
  air.ToMoleFraction();
  //   for (int k = ks; k <= ke; ++k)
    // for (int j = js; j <= je; ++j) {
  air.SetZero();
  air.w[iH2O] = xH2O;
  air.w[iNH3] = xNH3;
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;

  int i = is;
  for (; i <= ie; ++i) {
    // check relative humidity
    Real rh = get_relative_humidity(air, iNH3);
    air.w[iNH3] *= std::min(rh_max_nh3 / rh, 1.);

    AirParcelHelper::distribute_to_primitive(pmb, ks, js, i, air);

    pthermo->Extrapolate(&air, -dlnp, "dry");

    if (air.w[IDN] < Tmin) break;
  }

  // Replace adiabatic atmosphere with isothermal atmosphere if temperature
  // is too low
  pthermo->Extrapolate(&air, dlnp, "dry");
  for (; i <= ie; ++i) {
    pthermo->Extrapolate(&air, -dlnp, "isothermal");
    AirParcelHelper::distribute_to_primitive(pmb, ks, js, i, air);
  }
    // }
};
