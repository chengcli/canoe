//! \file nbody_inbound.cpp
//! \brief Implementations of inbound functions to particle

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/interpolation.h>  // interpn

// particle
#include "particle_data.hpp"
#include "particles.hpp"

void ParticleBase::SetVelocitiesFromHydro(Hydro const *phydro,
                                          Coordinates const *pcoord) {
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];

  v1.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(phydro->w), 4, IM1,
                          1);
  v2.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(phydro->w), 4, IM2,
                          1);
  v3.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(phydro->w), 4, IM3,
                          1);

  auto pm = pmy_block_->pmy_mesh;

  // loop over particles
  for (auto &q : pc) {
    loc[0] = q.x3;
    loc[1] = q.x2;
    loc[2] = q.x1;

    int f = pm->f2 + pm->f3;
    // interpolate Eulerian velocity to particle velocity
    interpn(&q.v1, loc + 2 - f, v1.data(), pcoord->GetCellCoords() + 2 - f,
            pcoord->GetDimensions() + 2 - f, 1 + f, 1);

    if (pm->f2) {
      interpn(&q.v2, loc + 2 - f, v2.data(), pcoord->GetCellCoords() + 2 - f,
              pcoord->GetDimensions() + 2 - f, 1 + f, 1);
    } else {
      q.v2 = 0.;
    }

    if (pm->f3) {
      interpn(&q.v3, loc + 2 - f, v3.data(), pcoord->GetCellCoords() + 2 - f,
              pcoord->GetDimensions() + 2 - f, 1 + f, 1);
    } else {
      q.v3 = 0.;
    }
  }
}
