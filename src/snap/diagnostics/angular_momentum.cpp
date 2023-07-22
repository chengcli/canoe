/** @file angular_momentum.cpp
 * @brief Implement total mass, moment of inertia, axial angular momentum output
 * (am)
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 18, 2021 08:07:22 PDT
 * @attention This output only works with spherical polar coordinate system
 * @bug No known bugs.
 */

// C/C++ header
#include <cstring>
#include <sstream>
#include <stdexcept>
// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Athena++ header
#include "../coordinates/coordinates.hpp"
#include "../math/core.h"
#include "diagnostics.hpp"

AngularMomentum::AngularMomentum(MeshBlock *pmb) : Diagnostics(pmb, "am") {
  type = "VECTORS";
  grid = "---";
  long_name =
      "mass,moment of inertia relative to a thin spherical shell,mean angular "
      "velocity";
  units = "kg,1,1/s";
  data.NewAthenaArray(3, 1, 1, 1);
  if (strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in AngularMomentum::AngularMomentum" << std::endl
        << "Diagnostics 'am' can only be used in spherical_polar geometry."
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

void AngularMomentum::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real mass_moi_am[3];
  mass_moi_am[0] = 0.;
  mass_moi_am[1] = 0.;
  mass_moi_am[2] = 0.;

  Real xmin = pmb->pmy_mesh->mesh_size.x1min;
  Real xmax = pmb->pmy_mesh->mesh_size.x1max;
  Real mean_radius = (xmin + xmax) / 2.;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        // mass
        mass_moi_am[0] += vol_(i) * w(IDN, k, j, i);
        // moment of inertia
        mass_moi_am[1] += vol_(i) * w(IDN, k, j, i) *
                          sqr(pcoord->x1v(i) * sin(pcoord->x2v(j)));
        // angular momentum
        mass_moi_am[2] += vol_(i) * w(IDN, k, j, i) * pcoord->x1v(i) *
                          sin(pcoord->x2v(j)) * w(IV3, k, j, i);
      }
    }

    // sum over all ranks
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, mass_moi_am, 3, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  data(0) = mass_moi_am[0];
  data(1) =
      mass_moi_am[1] / (2. / 3. * mass_moi_am[0] * mean_radius * mean_radius);
  if (mass_moi_am[1] != 0) data(2) = mass_moi_am[2] / mass_moi_am[1];
}
