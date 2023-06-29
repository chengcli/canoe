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

// canoe
#include <configure.hpp>
#include <impl.hpp>
#include <constants.hpp>
#include <variable.hpp>
#include <index_map.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics_helper.hpp>
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


Real grav, P0, T0, xHe, xCH4, Tmin, clat;
Real xH2S, xNa, xKCl, metallicity;

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
  auto pthermo = pimpl->pthermo;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->PotentialTemp(phydro->w.at(k,j,i), P0);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        // mse
        user_out_var(3,k,j,i) = pthermo->MoistStaticEnergy(phydro->w.at(k,j,i), grav*pcoord->x1v(i));
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3+n,k,j,i) = pthermo->RelativeHumidity(phydro->w.at(k,j,i), n);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("mesh", "ReferencePressure");
  T0 = pin->GetReal("problem", "T1bar");

  Tmin = pin->GetReal("problem", "Tmin");
  clat = pin->GetOrAddReal("problem", "clat", 0.);

  //xH2S = pin->GetReal("problem", "xH2S");
  //metallicity = pin->GetOrAddReal("problem", "metallicity", 0.);
  //xNa = pin->GetReal("problem", "xNa");
  //xNa *= pow(10., metallicity);
  //xKCl = pin->GetReal("problem", "xKCl");
  //xKCl *= pow(10., metallicity);
}

void update_gamma(Real* gamma, Real const q[]) {
  //std::cout << "I'm here" << std::endl;
  Real T = q[IDN], cp_h2, cp_he, cp_ch4;
  if (T < 300.)
    cp_h2 = Hydrogen::cp_norm(T);
  else
    cp_h2 = Hydrogen::cp_nist(T);
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

  Real cp_real = (1. - xHe - xCH4)*cp_h2 + xHe*cp_he + xCH4*cp_ch4;
  *gamma = cp_real/(cp_real - Constants::Rgas);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Application::Logger app("main");
  app->Log("ProblemGenerator: juno");

  // mesh limits
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real H0 = pimpl->GetPressureScaleHeight();

  // request temperature and pressure
  app->Log("request T", T0);
  app->Log("request P", P0);

  // thermodynamic constants
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma/(gamma - 1.)*Rd;

  // estimate surface temperature and pressure
  Real Ps = P0*exp(-x1min/H0);
  Real Ts = T0*pow(Ps/P0,Rd/cp);
  Real dlnp, **w1, *z1, *t1;
  dlnp = (x1max - x1min)/pmy_mesh->mesh_size.nx1/(2.*H0);
  size_t nx1 = (x1max - x1min)/(H0*dlnp);
  NewCArray(w1, nx1, Variable::Size);
  z1 = new Real [nx1];
  t1 = new Real [nx1];

  // index
  int iH2O = pindex->GetVaporId("H2O");
  int iNH3 = pindex->GetVaporId("NH3");

  // read in vapors
  w1[0][iH2O] = pin->GetReal("problem", "qH2O.gkg")/1.E3;
  w1[0][iNH3] = pin->GetReal("problem", "qNH3.gkg")/1.E3;
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + H0*dlnp;

  auto pthermo = pimpl->pthermo;

  // set up an adiabatic atmosphere
  int max_iter = 200, iter = 0;
  Real t0;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, -dlnp, nx1, Adiabat::dry, 1.);

    // Replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1-1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // Find T at p = p0
    for (int i = 0; i < nx1; ++i) {
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    t0 = interp1(0, t1, z1, nx1);

    Ts += T0 - t0;
    if (fabs(T0 - t0) < 0.01) break;

    app->Log("iteration #", iter);
    app->Log("T", t0);
  }

  if (iter > max_iter) {
    throw InvalidValueError("ProblemGenerator", "maximum iteration reached");
  }

  // change to log quantity
  for (int i = 0; i < nx1; ++i)
    for (int n = 0; n < Variable::Size; ++n) {
      if (n == IVX || n == IVY || n == IVZ)
        w1[i][n] = 0.;
      else
        w1[i][n] = log(w1[i][n]);
    }

  // interpolate to model grid
  for (int i = is; i <= ie; ++i) {
    Real buf[Variable::Size];

    // set gas concentration and velocity
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, Variable::Size);
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) {
          if (n == IVX || n == IVY || n == IVZ)
            phydro->w(n,k,j,i) = 0.;
          else {
            phydro->w(n,k,j,i) = exp(buf[n]);
          }
        }

    // set tracers, electron, Na and K
    int ielec = pindex->GetTracerId("e-");
    int iNa = pindex->GetTracerId("Na");
    int iK = pindex->GetTracerId("K");
    auto ptracer = pimpl->ptracer;

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real temp = pthermo->GetTemp(phydro->w.at(k,j,i));
        Real pH2S = xH2S*phydro->w(IPR,k,j,i);
        Real pNa = xNa*phydro->w(IPR,k,j,i);
        Real svp = sat_vapor_p_Na_H2S_Visscher(temp, pH2S);
        pNa = std::min(svp, pNa);

        ptracer->s(iNa,k,j,i) = pNa/(Constants::kBoltz*temp);
        ptracer->s(ielec,k,j,i) = saha_ionization_electron_density(temp, ptracer->s(iNa,k,j,i), 5.14);
      }
  }

  // if requested
  // replace adiabatic temperature profile with measurement
  if (pin->GetOrAddBoolean("problem", "use_fletcher16_cirs", false)) {
    // clat->glat, pa -> bar
    Real glat = jup_centric2graphic(clat);
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          if (phydro->w(IPR,k,j,i) > 1.E5) continue;
          Real temp_ad = pthermo->GetTemp(phydro->w.at(k,j,i));
          // pa -> bar
          Real temp_real = Jupiter::get_temp_fletcher16_cirs(glat, phydro->w(IPR,k,j,i)/1.E5);
          Real R = pthermo->RovRd(phydro->w.at(k,j,i))*Rd;
          phydro->w(IDN,k,j,i) = phydro->w(IPR,k,j,i)/(R*std::max(temp_real, temp_ad));
        }
  }

  // primitive to conserved conversion
  peos->PrimitiveToConserved(phydro->w, pfield->bcc,
    phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] t1;
}

int main(int argc, char **argv)  {
  auto app = Application::GetInstance();

  app->InstallMonitor("main", "main.out", "main.err");
  app->InstallMonitor("inversion", "main.out", "main.err");
  app->InstallMonitor("astro", "main.out", "main.err");
  app->InstallMonitor("snap", "main.out", "main.err");
  app->InstallMonitor("harp", "main.out", "main.err");

  Application::Start();

  Application::Logger log("main");

  // no MPI
  Globals::my_rank = 0;
  Globals::nranks  = 1;

  IOWrapper infile;
  infile.Open("juno.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput;
  pinput->LoadFromFile(infile);
  infile.Close();

  // set up mesh
  int restart = false;
  int mesh_only = false;

  Mesh *pmesh = new Mesh(pinput, mesh_only);

  // set up components
  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock *pmb = pmesh->my_blocks(b);
    pmb->pindex = std::make_shared<MeshBlock::IndexMap>(pmb, pinput);
    pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
  }

  // initialize mesh
  pmesh->Initialize(restart, pinput);

  // calculate radiative transfer
  for (int b = 0; b < pmesh->nblocal; ++b) {
    auto pmb = pmesh->my_blocks(b);
    auto phydro = pmb->phydro;
    auto prad = pmb->pimpl->prad;

    int is = pmb->is, js = pmb->js, ks = pmb->ks;
    int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

    AthenaArray<Real> tempa, qNH3a, qH2Oa;
    // read profile updates from input
    std::string file_tempa = pinput->GetOrAddString("problem", "file.tempa", "");
    if (file_tempa != "") {
      ReadDataTable(&tempa, file_tempa);
    }
    std::string file_nh3a = pinput->GetOrAddString("problem", "file.nh3a", "");
    if (file_nh3a != "") {
      ReadDataTable(&qNH3a, file_nh3a);
    }
    std::string file_h2oa = pinput->GetOrAddString("problem", "file.h2oa", "");
    if (file_h2oa != "") {
      ReadDataTable(&qH2Oa, file_h2oa);
    }

    for (int k = ks; k <= ke; ++k) {
      // run RT models
      log->Log("run microwave radiative transfer");
      for (int j = js-1; j <= je; ++j)
        prad->calculateRadiance(&prad->radiance, 0., k, j, is, ie+1);
    }

  }

  // make output
  Outputs *pouts;
  pouts = new Outputs(pmesh, pinput);
  pouts->MakeOutputs(pmesh, pinput);

  delete pinput;
  delete pmesh;
  delete pouts;

  Application::Destroy();
}
