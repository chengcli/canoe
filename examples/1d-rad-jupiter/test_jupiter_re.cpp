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

// utils
#include <utils/ndarrays.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/thermodynamics_helper.hpp>

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
  auto pthermo = pimpl->pthermo;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0, k, j, i) =
            pthermo->GetTemp(phydro->w.at(k, j, i));
        user_out_var(1, k, j, i) =
            pthermo->PotentialTemp(phydro->w.at(k, j, i), P0);
      }
}

/*void RadiationBand::addAbsorber(
        std::string name, std::string file,
                                ParameterInput *pin) {
  Real xHe = pin->GetOrAddReal("radiation", "xHe", 0.136);
  Real xH2 = 1. - xHe;

  std::stringstream msg;

  if (name == "H2-H2") {
    pabs->addAbsorber(XizH2H2CIA(this, 0, xH2))->loadCoefficient(file);
  } else if (name == "H2-He") {
    pabs->addAbsorber(XizH2HeCIA(this, 0, xH2, xHe))->loadCoefficient(file);
  } else if (name == "freedman_simple") {
    pabs->addAbsorber(FreedmanSimple(this, pin));
  } else if (name == "freedman_mean") {
    pabs->addAbsorber(FreedmanMean(this));
  } else if (strncmp(name.c_str(), "ck-", 3) == 0) {
    char str[80], aname[80];
    int bid;
    strcpy(str, name.c_str());
    replaceChar(str, '-', ' ');
    int err = sscanf(str, "ck %s %d", aname, &bid);
    if (err != EOF) {
      pabs->addAbsorber(CorrelatedKAbsorber(this, aname))
          ->loadCoefficient(file, bid);
    } else {
      msg << "### FATAL ERROR in RadiationBand::addAbsorber" << std::endl
          << "Incorrect format for absorber '" << name << "' ";
      ATHENA_ERROR(msg);
    }
  } else {
    msg << "### FATAL ERROR in RadiationBand::addAbsorber" << std::endl
        << "Unknown absorber: '" << name << "' ";
    ATHENA_ERROR(msg);
  }
}*/

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
  auto pthermo = pimpl->pthermo;
  Real Rd = pthermo->GetRd();
  Real cp = gamma / (gamma - 1.) * Rd;
  Real Ts = T0 - grav / cp * (x1min - Z0);
  Real Ps = P0 * pow(Ts / T0, cp / Rd);
  int max_iter = 200, iter = 0;

  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i) z1[i] = z1[i - 1] + dz;

  Real t0, p0;
  if (Globals::my_rank == 0)
    std::cout << "- request T = " << T0 << " P = " << P0 << " at Z = " << Z0
              << std::endl;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo,
                                 0.);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if
    // temperature is too low
    int ii = 0;
    for (; ii < nx1 - 1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR] / (w1[ii][IDN] * Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR] * exp(-grav * (z1[i] - z1[ii]) / (Rd * Tv));
      w1[i][IDN] = w1[i][IPR] / (Rd * Tv);
      for (int n = 1; n <= NVAPOR; ++n) w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z = Z0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    p0 = interp1(Z0, p1, z1, nx1);
    t0 = interp1(Z0, t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0 / p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0 / p0 - 1.) < 1.E-4)) break;
    if (Globals::my_rank == 0)
      std::cout << "- iteration #" << iter << ": "
                << "T = " << t0 << " P = " << p0 << std::endl;
  }

  if (iter > max_iter) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator" << std::endl
        << "maximum iteration reached." << std::endl
        << "T0 = " << t0 << std::endl
        << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  int kl = block_size.nx3 == 1 ? ks : ks - NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke + NGHOST;
  int jl = block_size.nx2 == 1 ? js : js - NGHOST;
  int ju = block_size.nx2 == 1 ? je : je + NGHOST;

  unsigned int seed = time(NULL);

  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;

    // set gas concentration
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) phydro->w(n, k, j, i) = buf[n];

    // add noise
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        phydro->w(IVX, k, j, i) = 0.01 * (1. * rand_r(&seed) / RAND_MAX - 0.5);
  }

  // set spectral properties
  auto prad = pimpl->prad;
  for (int b = 0; b < prad->GetNumBands(); ++b) {
    auto p = prad->GetBand(b);

    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        p->SetSpectralProperties(k, j, is, ie);
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}

int main(int argc, char **argv) {
#ifndef HYDROSTATIC
  static_assert(false, "This problem requires turning on hydrostatic option");
#endif
  std::stringstream msg;

  // no MPI
  Globals::my_rank = 0;
  Globals::nranks = 1;

  IOWrapper infile;
  infile.Open("saturn_atmosphere.inp", IOWrapper::FileMode::read);

  ParameterInput *pinput = new ParameterInput;
  pinput->LoadFromFile(infile);
  infile.Close();

  // set up mesh
  int restart = false;
  int mesh_only = false;

  auto app = Application::GetInstance();

  auto logger = app->GetMonitor("main");

  Mesh *pmesh = new Mesh(pinput, mesh_only);

  // set up components
  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock *pmb = pmesh->my_blocks(b);
    pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
  }

  // initialize mesh
  pmesh->Initialize(restart, pinput);
  auto pindex = pmesh->my_blocks(0)->pindex;
  auto psnap = pmesh->my_blocks(0)->pimpl;

  // mesh limits
  Real x1min = pmesh->mesh_size.x1min;
  Real x1max = pmesh->mesh_size.x1max;
  Real H0 = psnap->GetPressureScaleHeight();

  // request temperature and pressure
  logger->Log("request T = " + std::to_string(T0));
  logger->Log("request P = " + std::to_string(P0));

  // thermodynamic constants
  Real gamma = pinput->GetReal("hydro", "gamma");
  Real Rd = pinput->GetReal("thermodynamics", "Rd");
  Real cp = gamma / (gamma - 1.) * Rd;
  Real Tmin = pinput->GetReal("problem", "Tmin");

  // estimate surface temperature and pressure
  Real Ps = P0 * exp(-x1min / H0);
  Real Ts = T0 * pow(Ps / P0, Rd / cp);
  Real dlnp, **w1, *z1, *t1;
  dlnp = (x1max - x1min) / pmesh->mesh_size.nx1 / (2. * H0);
  size_t nx1 = static_cast<size_t>((x1max - x1min) / (H0 * dlnp));
  NewCArray(w1, nx1, NHYDRO + 2 * NVAPOR);
  z1 = new Real[nx1];
  t1 = new Real[nx1];

  // read in vapors
  size_t iH2O = pindex->GetVaporId("H2O");
  size_t iNH3 = pindex->GetVaporId("NH3");

  w1[0][iH2O] = pinput->GetReal("problem", "qH2O");
  w1[0][iNH3] = pinput->GetReal("problem", "qNH3");
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i) z1[i] = z1[i - 1] + H0 * dlnp;

  // set up variables
  for (int b = 0; b < pmesh->nblocal; ++b) {
    MeshBlock *pmb = pmesh->my_blocks(b);
    auto phydro = pmb->phydro;
    auto pfield = pmb->pfield;
    auto pcoord = pmb->pcoord;
    auto peos = pmb->peos;
    auto pthermo = pmb->pimpl->pthermo;

    // set up an adiabatic atmosphere
    int max_iter = 200, iter = 0;
    Real t0;
    while (iter++ < max_iter) {
      pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, -dlnp, nx1,
                                   Adiabat::pseudo, 1.);

      // Replace adiabatic atmosphere with isothermal atmosphere if temperature
      // is too low
      int ii = 0;
      for (; ii < nx1 - 1; ++ii)
        if (pthermo->GetTemp(w1[ii]) < Tmin) break;
      Real Tv = w1[ii][IPR] / (w1[ii][IDN] * Rd);
      for (int i = ii; i < nx1; ++i) {
        w1[i][IDN] = w1[i][IPR] / (Rd * Tv);
        for (int n = 1; n <= NVAPOR; ++n) w1[i][n] = w1[ii][n];
      }

      // Find T at p = p0
      for (int i = 0; i < nx1; ++i) {
        t1[i] = pthermo->GetTemp(w1[i]);
      }
      t0 = interp1(0, t1, z1, nx1);

      Ts += T0 - t0;
      if (fabs(T0 - t0) < 0.01) break;

      logger->Log("Iteration # = " + std::to_string(iter));
      logger->Log("T = " + std::to_string(t0));
    }

    if (iter > max_iter) {
      msg << "### FATAL ERROR" << std::endl
          << "maximum iteration reached." << std::endl
          << "T0 = " << t0;
      ATHENA_ERROR(msg);
    }

    // Change to log quantity
    for (int i = 0; i < nx1; ++i)
      for (int n = 0; n < NHYDRO + 2 * NVAPOR; ++n) {
        if (n == IVX || n == IVY || n == IVZ)
          w1[i][n] = 0.;
        else
          w1[i][n] = log(w1[i][n]);
      }

    // interpolate to model grid
    int is = pmb->is, js = pmb->js, ks = pmb->ks;
    int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real buf[NHYDRO + 2 * NVAPOR];

          // set gas concentration and velocity
          interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO + 2 * NVAPOR);
          for (int n = 0; n < NHYDRO; ++n)
            for (int k = ks; k <= ke; ++k)
              for (int j = js; j <= je; ++j) {
                if (n == IVX || n == IVY || n == IVZ) {
                  phydro->w(n, k, j, i) = 0.;
                } else {
                  phydro->w(n, k, j, i) = exp(buf[n]);
                }
              }
        }

    // primitive to conserved conversion
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pmb->pcoord,
                               is, ie, js, je, ks, ke);
  }

  // make output
  Outputs *pouts;
  pouts = new Outputs(pmesh, pinput);
  pouts->MakeOutputs(pmesh, pinput);

  FreeCArray(w1);
  delete[] z1;
  delete[] t1;

  delete pinput;
  delete pmesh;
  delete pouts;
}
