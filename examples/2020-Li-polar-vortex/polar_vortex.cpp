// changes: no time loop, exp v2, vphi = vmax (positive #), coriolis not
// approximated, damping model--Model #2--ingersoll's equation n=2
//! \file polar_vortex.cpp
//  \brief jupiter polar vortex model

// C/C++
#include <cmath>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

// boost
#include <boost/math/special_functions/gamma.hpp>

// Athena++
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/core.h>  // sqr

// utils
#include <utils/fileio.hpp>  // ReadDataTable

// global parameters
Real phi0, phi1, lambda, alpha, vphi, vrad, vis, polarity, skewness, omega,
    radius, sponge_lat, sponge_tau, max_cyclone, max_anticyclone, b;

Real tseed, interval;

std::default_random_engine gen;
std::uniform_real_distribution<double> uniform(0.0, 1.0);

void FindLatlon(Real *lat, Real *lon, Real x2, Real x1) {
  Real dist = sqrt(x1 * x1 + x2 * x2);
  *lat = M_PI / 2. - dist / radius;
  *lon = asin(x1 / dist);
  if (x2 > 0 && x1 > 0) *lon = M_PI - *lon;
  if (x2 > 0 && x1 < 0) *lon = -M_PI - *lon;
}

void FindXY(Real *x1, Real *x2, Real lat, Real lon) {
  Real dist = radius * (M_PI / 2. - lat);
  *x1 = dist * sin(lon);
  *x2 = -dist * cos(lon);
}

Real Phib(Coordinates *pcoord, int j,
          int i)  // phi1 determines topography at the pole
{
  Real dist = sqrt(sqr(pcoord->x1v(i)) + sqr(pcoord->x2v(j)));
  return phi1 * exp(-pow(dist / lambda, alpha) / alpha);
}

void PolarVortexForcing(MeshBlock *pmb, const Real time, const Real dt,
                        AthenaArray<Real> const &w, AthenaArray<Real> const &r,
                        AthenaArray<Real> const &bcc, AthenaArray<Real> &u,
                        AthenaArray<Real> &s) {
  for (int j = pmb->js; j <= pmb->je; ++j)
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real x1 = pmb->pcoord->x1v(i);
      Real x2 = pmb->pcoord->x2v(j);
      Real dist = sqrt(x1 * x1 + x2 * x2);

      // coriolis force
      Real f = 2. * omega * cos(dist / radius);  // not approximated
      u(IM1, j, i) += dt * f * w(IDN, j, i) * w(IVY, j, i);
      u(IM2, j, i) += -dt * f * w(IDN, j, i) * w(IVX, j, i);

      // topography
      Real phib = Phib(pmb->pcoord, j, i);
      u(IM1, j, i) +=
          dt * w(IDN, j, i) * phib * x1 * pow(dist / lambda, alpha) / sqr(dist);
      u(IM2, j, i) +=
          dt * w(IDN, j, i) * phib * x2 * pow(dist / lambda, alpha) / sqr(dist);

      Real s = (90. - dist / radius / M_PI * 180.) / sponge_lat;
      if (s < 1) {  // sponge layer
        // u(IM1,j,i) -= dt*u(IM1,j,i)/pow(sponge_tau,s);
        // u(IM2,j,i) -= dt*u(IM1,j,i)/pow(sponge_tau,s);
        // Real lat, lon;
        // FindLatlon(&lat, &lon, x2, x1);
        u(IM1, j, i) -=
            dt * w(IDN, j, i) * w(IVX, j, i) * pow((1. - s), 2) / sponge_tau;
        u(IM2, j, i) -=
            dt * w(IDN, j, i) * w(IVY, j, i) * pow((1. - s), 2) / sponge_tau;
        // u(ITR,j,i) -= dt*(u(ITR,j,i) - lat/M_PI*180.)/pow(sponge_tau,s);
      } else {  // viscosity
        Real dx = pmb->pcoord->dx1f(i);
        Real dy = pmb->pcoord->dx2f(j);
        u(IM1, j, i) +=
            vis * dt / (dx * dy) * w(IDN, j, i) *
            (w(IM1, j + 1, i) + w(IM1, j - 1, i) + w(IM1, j, i + 1) +
             w(IM1, j, i - 1) - 4 * w(IM1, j, i));
        u(IM2, j, i) +=
            vis * dt / (dx * dy) * w(IDN, j, i) *
            (w(IM2, j + 1, i) + w(IM2, j - 1, i) + w(IM2, j, i + 1) +
             w(IM2, j, i - 1) - 4 * w(IM2, j, i));
      }
    }
}
/*
void MeshBlock::UserWorkInLoop()        //called at the end of every time step
{
    Real time = pmy_mesh->time;
    // check if need to seed a new vortex
    if ((time > 100) && (time < 400)) {
        std::cout << "UWIL time: " << time << "\n";
        Real x1, x2, range;
        //Real min_lat = sponge_lat/180.*M_PI + 1.5*vrad/radius;
        //Real vpol = uniform(gen) > polarity ? 1 : -1;       //? makes an if
else statement, uniform(gen) makes random number Real vpol = 1;
        //if (vpol > 0) // cyclone        //basically a 50/50 chance of
generating a cyclone or anticyclone? depends on dist
        //range = std::max(0., sin(max_cyclone/180.*M_PI) - sin(min_lat));
        //else { // anticyclone
            //range = std::max(0., sin(max_anticyclone/180.*M_PI) -
sin(min_lat));
            //vpol *= skewness;
        //}
        //Real vlat = asin(sin(min_lat) + uniform(gen)*range);    //generates a
random lat and lon in the computed range Real vlat = uniform(gen) * 15. + 75.;
        Real vlon = uniform(gen)*2.*M_PI;
        std::cout << "lat and lon: " << vlat << " " << vlon << "\n";
        FindXY(&x1, &x2, vlat, vlon);       //then finds corresponding xy
        for (int j = js - NGHOST; j <= je + NGHOST; ++j)
        for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
            Real lat, lon, vel;
            FindLatlon(&lat, &lon, pcoord->x2v(j), pcoord->x1v(i));
            Real dist = sqrt(sqr(x1 - pcoord->x1v(i)) + sqr(x2 -
pcoord->x2v(j))); Real phi = vpol*(phi0/2)*exp(-0.5*sqr(dist)/sqr(500000));
//sets pressure field Real fcor = 2.*omega*sin(lat); if (vpol > 0){  // cyclone
            vel = -dist*fcor/2. + dist*fcor/2.*sqrt(1.
- 4.*phi/(sqr(fcor*500000)));     //cyclostrophic vel. magnitude } else { //
anticyclone vel = -phi*dist/(fcor*vrad*vrad); //why is this different?
            }
            //Real vel = -dist*fcor/2. + dist*fcor/2.*sqrt(1.
- 4.*phi/(sqr(fcor*vrad))); phydro->w(IDN,j,i) += phi; phydro->w(IVX,j,i) +=
vel*(x2 - pcoord->x2v(j))/dist; phydro->w(IVY,j,i) += -vel*(x1 -
pcoord->x1v(i))/dist;
        }
        peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord,
        is - NGHOST, ie + NGHOST, js - NGHOST, je + NGHOST, ks, ke);
        tseed = time - interval*log(uniform(gen));          //sets a random time
until the next one will be generated. but in interval = 9E9 ???
    }
}
*/
void Mesh::InitUserMeshData(
    ParameterInput *pin)  // called at beginning of simulation to initialize
                          // global parameters
{
  // forcing parameters
  phi0 = pin->GetReal("problem", "phi0");
  phi1 = pin->GetReal("problem", "phi1");
  lambda = pin->GetReal("problem", "lambda");
  alpha = pin->GetOrAddReal("problem", "alpha", 1.);
  radius = pin->GetReal("problem", "radius");
  sponge_lat = pin->GetReal("problem", "sponge_lat");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  omega = pin->GetReal("problem", "omega");
  interval = pin->GetReal("problem", "interval");
  vrad = pin->GetReal("problem", "vrad");
  vphi = pin->GetReal("problem", "vphi");
  vis = pin->GetOrAddReal("problem", "vis", 0.);
  max_cyclone = pin->GetOrAddReal("problem", "max_cyclone", 90.);
  max_anticyclone = pin->GetOrAddReal("problem", "max_anticyclone", 90.);
  polarity = pin->GetOrAddReal("problem", "polarity", 0.);
  skewness = pin->GetOrAddReal("problem", "skewness", 1.);
  b = 2.;

  // next seeding time
  tseed = time;
  // tseed = time - interval*log(uniform(gen));

  // forcing function
  EnrollUserExplicitSourceFunction(PolarVortexForcing);
}

void MeshBlock::InitUserMeshBlockData(
    ParameterInput *pin)  // allocating and initializing local variables
                          // belonging to each MeshBlock
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "vort");
  SetUserOutputVariableName(1, "pv");
}

void MeshBlock::UserWorkBeforeOutput(
    ParameterInput *pin)  // called right before each output file is produced
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dx = pcoord->dx1f(i);
        Real dy = pcoord->dx2f(j);
        // user_out_var(0,k,j,i) = (phydro->vf(1,k,j,i+1) -
        // phydro->vf(1,k,j,i))/dx
        //                       - (phydro->vf(3,k,j+1,i) -
        //                       phydro->vf(3,k,j,i))/dy;
        user_out_var(0, k, j, i) =
            (phydro->w(IVY, k, j, i + 1) - phydro->w(IVY, k, j, i - 1)) /
                (2. * dx) -
            (phydro->w(IVX, k, j + 1, i) - phydro->w(IVX, k, j - 1, i)) /
                (2. * dy);
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 2. * omega * cos(dist / radius);  // not approx
        // Real f = 2.*omega;
        user_out_var(1, k, j, i) =
            (f + user_out_var(0, k, j, i)) / phydro->w(IDN, k, j, i);
      }
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(
    ParameterInput *pin)  // sets initial conditions
{
  AthenaArray<Real> vpos;
  ReadDataTable(&vpos, pin->GetString("problem", "vfile"));
  int vnum = vpos.GetDim2();  // input file allows positions
  // int vnum = 6;
  for (int i = 0; i < vnum; ++i) {
    vpos(i, 0) *= M_PI / 180.;  // convert to radians
    vpos(i, 1) *= M_PI / 180.;
  }

  // setup initial height and velocity field
  for (int j = js - NGHOST; j <= je + NGHOST; ++j)
    for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
      Real lat, lon;
      FindLatlon(&lat, &lon, pcoord->x2v(j), pcoord->x1v(i));

      // background
      phydro->w(IDN, j, i) = phi0 - Phib(pcoord, j, i);  // set background
      // phydro->w(IDN,j,i) = phi0;

      // add contribution from vortices
      for (int n = 0; n < vnum; ++n) {
        Real x1, x2, vel;
        Real fcor = 2. * omega * sin(lat);  // no approx
        // Real fcor = 2.*omega;
        FindXY(&x1, &x2, vpos(n, 0), vpos(n, 1));
        Real dist = sqrt(sqr(x1 - pcoord->x1v(i)) + sqr(x2 - pcoord->x2v(j)));
        Real vsize = vrad * vpos(n, 2);
        // Real phi = vpos(n,3)*vphi*exp(-0.5*sqr(dist)/sqr(vsize)); //original
        // Real phi = vpos(n,3)*(-vphi*fcor*(dist + vsize)*exp(1 - (dist/vsize))
        // - (sqr(vphi)/(4*vsize))*(2*dist + vsize)*exp(2-((2*dist)/vsize)));
        // //b=1 two term force balance

        Real gam = boost::math::tgamma(2. / b, (1. / b) * pow(dist / vsize, b));

        Real phi = vpos(n, 3) * -fcor * vphi * vsize * exp(1. / b) *
                   pow(b, 2. / b - 1.) * gam;
        // vpos(n,3) is the polarity and thus determines the sign!

        if (vpos(n, 3) > 0) {  // cyclone
                               // vel = -dist*fcor/2. + dist*fcor/2.*sqrt(1.
          // - 4.*phi/(sqr(fcor*vsize)));        //original vel =
          // (vphi/vsize)*dist*exp(1-(dist/vsize));  //mine v2
          vel = (vphi / vsize) * dist *
                exp((1. / b) * (1 - pow((dist / vsize), b)));
        } else {  // anticyclone
          vel = -(vphi / vsize) * dist *
                exp((1. / b) * (1 - pow((dist / vsize), b)));
        }
        phydro->w(IDN, j, i) += phi;
        phydro->w(IVX, j, i) += vel * (x2 - pcoord->x2v(j)) / dist;
        phydro->w(IVY, j, i) += -vel * (x1 - pcoord->x1v(i)) / dist;
      }
    }

  /* add tracer
   for (int j = js; j <= je; ++j)
   for (int i = is; i <= ie; ++i) {
   Real lat, lon;
   FindLatlon(&lat, &lon, pcoord->x2v(j), pcoord->x1v(i));
   phydro->w(ITR,j,i) = lat/M_PI*180.;
   }*/

  // convert to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie,
                             js, je, ks, ke);
}
