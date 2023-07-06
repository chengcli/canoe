// C/C++
#include <vector>
#include <cstdio>
#include <memory>
#include <fstream>
#include <algorithm>

// athena
#include <athena/parameter_input.hpp>
#include <athena/eos/eos.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/field/field.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/stride_iterator.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>
#include <constants.hpp>
#include <variable.hpp>
#include <index_map.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/molecules.hpp>
#include <snap/thermodynamics/vapors/sodium_vapors.hpp>

// utils
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>
#include <utils/fileio.hpp>

// astro
#include <astro/Jupiter/jup_fletcher16_cirs.hpp>
#include <astro/planets.hpp>

// harp
#include <harp/radiation_band.hpp>

// opacity
#include <opacity/Giants/microwave/mwr_absorbers.hpp>

// inversion
#include <inversion/profile_inversion.hpp>

Real grav, P0, T0, Tmin, clat;
Real xHe, xCH4, xH2S, xNa, xKCl, metallicity;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4+NVAPOR);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3+n, name.c_str());
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  auto pthermo = Thermodynamics::GetInstance();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(this,k,j,i);
        user_out_var(1,k,j,i) = pthermo->PotentialTemp(this,P0,k,j,i);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(this,k,j,i);
        // mse
        user_out_var(3,k,j,i) = pthermo->MoistStaticEnergy(this,grav*pcoord->x1v(i),k,j,i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3+n,k,j,i) = pthermo->RelativeHumidity(this,n,k,j,i);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("mesh", "ReferencePressure");
  T0 = pin->GetReal("problem", "T1bar");

  Tmin = pin->GetReal("problem", "Tmin");
  clat = pin->GetOrAddReal("problem", "clat", 0.);

  xH2S = pin->GetReal("problem", "xH2S");

  metallicity = pin->GetOrAddReal("problem", "metallicity", 0.);

  xNa = pin->GetReal("problem", "xNa");
  xNa *= pow(10., metallicity);

  xKCl = pin->GetReal("problem", "xKCl");
  xKCl *= pow(10., metallicity);

  xHe = pin->GetReal("problem", "xKCl");

  xCH4 = pin->GetReal("problem", "xKCl");
}

Real Thermodynamics::GetGammad(Variable const& qfrac) const {
  //std::cout << "I'm here" << std::endl;
  Real T = qfrac.w[IDN], cp_h2, cp_he, cp_ch4;
  if (T < 300.) {
    cp_h2 = Hydrogen::cp_norm(T);
  } else {
    cp_h2 = Hydrogen::cp_nist(T);
  }
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

  Real cp_real = (1. - xHe - xCH4)*cp_h2 + xHe*cp_he + xCH4*cp_ch4;
  return cp_real/(cp_real - Constants::Rgas);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Application::Logger app("main");
  app->Log("ProblemGenerator: juno");

  auto pthermo = Thermodynamics::GetInstance();

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real H0 = pimpl->GetPressureScaleHeight();

  // request temperature and pressure
  app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;

  // index
  int iH2O = pindex->GetVaporId("H2O");
  int iNH3 = pindex->GetVaporId("NH3");

  app->Log("index of H2O", iH2O);
  app->Log("index of NH3", iNH3);

  // set up an adiabatic atmosphere
  int max_iter = 200, iter = 0;
  Real dlnp = pcoord->dx1f(is)/H0;

  Variable var;
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
  var.SetZero();

  var.w[iH2O] = xH2O;
  var.w[iNH3] = xNH3;
  var.w[IPR] = Ps;
  var.w[IDN] = Ts;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      int i = is;
      for (; i <= ie; ++i) {
        pimpl->DistributePrimitive(var, k, j, i);
        pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::DryAdiabat);
        if (var.w[IDN] < Tmin) break;
      }

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
      for (; i <= ie; ++i) {
        pimpl->DistributePrimitive(var, k, j, i);
        pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::Isothermal);
      }
    }

  // set tracers, electron and Na
  int ielec = pindex->GetTracerId("e-");
  int iNa = pindex->GetTracerId("Na");
  auto ptracer = pimpl->ptracer;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real temp = pthermo->GetTemp(this,k,j,i);
        Real pH2S = xH2S*phydro->w(IPR,k,j,i);
        Real pNa = xNa*phydro->w(IPR,k,j,i);
        Real svp = sat_vapor_p_Na_H2S_Visscher(temp, pH2S);
        pNa = std::min(svp, pNa);

        ptracer->u(iNa,k,j,i) = pNa/(Constants::kBoltz*temp);
        ptracer->u(ielec,k,j,i) = saha_ionization_electron_density(temp, ptracer->u(iNa,k,j,i), 5.14);
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

        pimpl->GatherPrimitive(&var, k, j, ibegin);

        for (int i = ibegin; i < iend; ++i) {
          pthermo->Extrapolate(&var, -dlnp, Thermodynamics::Method::DryAdiabat,
              0., adlnTdlnP);
          pimpl->DistributePrimitive(var, k, j, i+1);
        }
      }
  }

  // replace adiabatic temperature profile with measurement
  if (pin->GetOrAddBoolean("problem", "use_fletcher16_cirs", false)) {
    // clat->glat, pa -> bar
    Real glat = jup_centric2graphic(clat);
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          if (phydro->w(IPR,k,j,i) > 1.E5) continue;
          Real temp_ad = pthermo->GetTemp(this,k,j,i);
          // pa -> bar
          Real temp_real = Jupiter::get_temp_fletcher16_cirs(glat, phydro->w(IPR,k,j,i)/1.E5);
          Real R = pthermo->RovRd(this,k,j,i)*Rd;
          phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)/(R*std::max(temp_real, temp_ad));
        }
  }

  // primitive to conserved conversion
  peos->PrimitiveToConserved(phydro->w, pfield->bcc,
    phydro->u, pcoord, is, ie, js, je, ks, ke);

  // Microwave radiative transfer needs temperatures at cell interfaces, which are
  // interpolated from cell centered hydrodynamic variables.
  // Normally, the boundary conditions are taken care of internally.
  // But, since we call radiative tranfer directly in pgen, we would need to update the
  // boundary conditions manually. The following lines of code updates the boundary
  // conditions.
  phydro->hbvar.SwapHydroQuantity(phydro->w, HydroBoundaryQuantity::prim);
  pbval->ApplyPhysicalBoundaries(0., 0., pbval->bvars_main_int);

  // calculate radiative transfer
  auto prad = pimpl->prad;

  for (int k = ks; k <= ke; ++k) {
    // run RT models
    app->Log("run microwave radiative transfer");
    for (int j = js; j <= je; ++j)
      prad->CalRadiance(0., k, j, is, ie+1);
  }
}
