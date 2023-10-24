// C/C++ header
#include <ctime>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../math/interpolation.h"
#include "../utils/utils.hpp" // replaceChar
#include "../thermodynamics/thermodynamics.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/hydrogen_cia.hpp"
#include "../radiation/freedman_mean.hpp"
#include "../radiation/freedman_simple.hpp"
#include "../radiation/correlatedk_absorber.hpp"
#include "../radiation/water_cloud.hpp"
#include "../physics/physics.hpp"
#include "../particles/particles.hpp"
#include "../communicator/communicator.hpp"
#include "../reconstruct/interpolation.hpp"
#include "../math/core.h"

// global parameters
enum {ic1 = 1}; //ic1 = 1, ic1c = NHYDRO+ic1-1, ic1p = NHYDRO+NVAPOR+ic1-1; 
Real grav, P0, T0, Z0, Tmin;
Real qvapor1, qvapor2, qRelaxT;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), P0, pthermo);
      }
}

void RadiationBand::addAbsorber(std::string name, std::string file, ParameterInput *pin)
{
  Real xHe = pin->GetOrAddReal("radiation", "xHe", 0.136);
  Real xH2 = 1. - xHe;

  std::stringstream msg;

  if (name == "H2O_l") {
    pabs->addAbsorber(FuWaterLiquidCloud(this, 0)); 
  } else if (name == "simplecloud") {
    pabs->addAbsorber(SimpleCloud(this, 0, pin));
  } else if (name == "H2-H2") {
    pabs->addAbsorber(XizH2H2CIA(this, 0, xH2))
        ->loadCoefficient(file);
  } else if (name == "H2-He") {
    pabs->addAbsorber(XizH2HeCIA(this, 0, xH2, xHe))
        ->loadCoefficient(file);
  } else if (name == "freedman_simple") {
    pabs->addAbsorber(FreedmanSimple(this, pin));
  } else if (name == "freedman_mean") {
    pabs->addAbsorber(FreedmanMean(this, pin));
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
      msg << "### FATAL ERROR in RadiationBand::addAbsorber"
          << std::endl << "Incorrect format for absorber '" << name <<"' ";
      ATHENA_ERROR(msg);
    }
  } else {
    msg << "### FATAL ERROR in RadiationBand::addAbsorber"
        << std::endl << "Unknown absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k){
    for (int j = js; j <= je; ++j){

      NeighborBlock const *pbot = pmb->pcomm->findBotNeighbor();
      if (pbot == nullptr) {
        // relax the vapor mixing ratio of the bottom cell to the initial value
        u(ic1,k,j,is) += w(IDN,k,j,is)*(qvapor1-w(ic1,k,j,is))*dt*qRelaxT;
      // assuming ocean below, so relax the vapor to the SVP

      }

//      for (int i = is; i <= ie; ++i) {
        // damp kinetic energy
//        u(IM1,k,j,i) -= w(IDN,k,j,i)*w(IM1,k,j,i)/5.;
//        u(IM2,k,j,i) -= w(IDN,k,j,i)*w(IM2,k,j,i)/5.;
//        u(IM3,k,j,i) -= w(IDN,k,j,i)*w(IM3,k,j,i)/5.;
//      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Z0 = pin->GetOrAddReal("problem", "Z0", 0.);
  Tmin = pin->GetReal("problem", "min_tem");
  qRelaxT = pin->GetReal("problem", "qRelaxT");
  qvapor1 = pin->GetReal("problem", "qvapor1");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  pdebug->Enter("ProblemGenerator: dry_rce");
  Real gamma = pin->GetReal("hydro", "gamma");

  // construct a 1D adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  Real **w1, *z1, *p1, *t1;
  Real dz = (x1max - x1min)/pmy_mesh->mesh_size.nx1/2.;
  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
//  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*(x1min - Z0);
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  for (int n = 1; n <= NVAPOR; ++n) { 
    Real qv = pin->GetReal("problem", "qvapor" + std::to_string(n));
    w1[0][n] = qv;
  }

  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  if (Globals::my_rank == 0)
    std::cout << "- request T = " << T0 << " P = " << P0 << " at Z = " << Z0 << std::endl;
  while (iter++ < max_iter) {
    //pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo, 0.);
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo, 1.E-3);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1-1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR]*exp(-grav*(z1[i] - z1[ii])/(Rd*Tv));
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z = Z0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    p0 = interp1(Z0, p1, z1, nx1);
    t0 = interp1(Z0, t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
    if (Globals::my_rank == 0)
      std::cout << "- iteration #" << iter << ": " << "T = " << t0 << " P = " << p0 << std::endl;
  }

  if (iter > max_iter) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
  int jl = block_size.nx2 == 1 ? js : js-NGHOST;
  int ju = block_size.nx2 == 1 ? je : je+NGHOST;
  srand(Globals::my_rank + time(0));
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
//    Real buf[NHYDRO];
//    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;

    // set gas concentration
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          phydro->w(n,k,j,i) = buf[n];

    // add noise
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
      phydro->w(IV1,k,j,i) = 0.01*(1.*rand()/RAND_MAX - 0.5);
  }

  // set spectral properties
  RadiationBand *p = prad->pband;
  while (p != NULL) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        p->setSpectralProperties(phydro->w, k, j, is, ie);
    p = p->next;
  }

  //pphy->Initialize(phydro->w);
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
  pdebug->Leave();
}
