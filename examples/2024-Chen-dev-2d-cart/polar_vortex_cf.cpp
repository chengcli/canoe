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

// global parameters
Real phi0, dphi, Li, vis, omega, sponge_lat, sponge_tau, max_freq, change_tau, radius;
int M;
AthenaArray<Real> As, Bs, Cs;

std::default_random_engine gen;
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);

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

      Real s = (90. - dist / radius / M_PI * 180.) / sponge_lat;
      if (s < 1) {
        u(IM1, j, i) -=
            dt * w(IDN, j, i) * w(IVX, j, i) * pow((1. - s), 2) / sponge_tau;
        u(IM2, j, i) -=
            dt * w(IDN, j, i) * w(IVY, j, i) * pow((1. - s), 2) / sponge_tau;
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

      // Forcing
      Real k = 2.0 * M_PI / Li;
      for (int m = 0; m < M; ++m) {
        Real phi = As(m) * dphi * cos(k * (x1 * cos(Cs(m)) + x2 * sin(Cs(m))) + Bs(m) * time);
        u(IDN, j, i) += dt * phi;
      }

    }
  // Evolve As, Bs, Cs
  for (int m = 0; m < M; ++m) {
    As(m) += dt * normal(gen) / change_tau;
    Bs(m) += dt * normal(gen) * max_freq / change_tau;
    Cs(m) += dt * normal(gen) * 2.0 * M_PI / change_tau;
  }
}


void Mesh::InitUserMeshData(
    ParameterInput *pin)  // called at beginning of simulation to initialize
                          // global parameters
{
  // forcing parameters
  phi0 = pin->GetReal("problem", "phi0");
  dphi = pin->GetReal("problem", "dphi");
  Li = pin->GetReal("problem", "Li");
  M = pin->GetInteger("problem", "M");
  sponge_lat = pin->GetReal("problem", "sponge_lat");
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  omega = pin->GetReal("problem", "omega");
  vis = pin->GetOrAddReal("problem", "vis", 0.);
  max_freq = pin->GetReal("problem", "max_freq");
  change_tau = pin->GetReal("problem", "change_tau");
  radius = pin->GetReal("problem", "radius");

  As.NewAthenaArray(M);
  Bs.NewAthenaArray(M);
  Cs.NewAthenaArray(M);
  for (int m = 0; m < M; ++m) {
    As(m) = 1.0;
    Bs(m) = uniform(gen) * max_freq;
    Cs(m) = uniform(gen) * 2.0 * M_PI;
  }

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
        user_out_var(0, k, j, i) =
            (phydro->w(IVY, k, j, i + 1) - phydro->w(IVY, k, j, i - 1)) /
                (2. * dx) -
            (phydro->w(IVX, k, j + 1, i) - phydro->w(IVX, k, j - 1, i)) /
                (2. * dy);
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real dist = sqrt(x1 * x1 + x2 * x2);
        Real f = 2. * omega * cos(dist / radius);

        user_out_var(1, k, j, i) =
            (f + user_out_var(0, k, j, i)) / phydro->w(IDN, k, j, i);
      }
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(
    ParameterInput *pin)  // sets initial conditions
{
  // setup initial height and velocity field
  for (int j = js - NGHOST; j <= je + NGHOST; ++j)
    for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
      Real lat, lon;
      FindLatlon(&lat, &lon, pcoord->x2v(j), pcoord->x1v(i));

      // background
      phydro->w(IDN, j, i) = phi0;  // set background

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
