// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>

// nbody
#include <nbody/particles/particles.hpp>

// canoe
#include "diagnostics.hpp"

Tendency::Tendency(MeshBlock *pmb) : Diagnostics(pmb, "tendency") {
  type = "VECTORS";
  grid = "--C";
  varname = "rho_t,";
  long_name = "horizontally averaged gas density tendency,";
  units = "kg/(m^3.s),";

  wh_.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
  for (int n = 1; n <= NVAPOR; ++n) {
    varname += "vapor" + std::to_string(n) + "_t,";
    units += "kg/(kg.s),";
    long_name += "horizontally averaged vapor mass mixing ratio tendency,";
  }
  for (int n = 0; n < 3; ++n) {
    units += "m/s^2,";
    varname += "vel" + std::to_string(n + 1) + "_t,";
    long_name += "horizontally averaged velocity tendency,";
  }
  units += "pa/s";
  varname += "press_t";
  long_name += "horizontally averaged pressure tendency,";

  int npart = 0;
  Particles *pp = pmb->ppart;

  while (pp != nullptr) {
    npart += pp->u.GetDim4();
    for (int n = 0; n < pp->u.GetDim4(); ++n) {
      units += ",kg/(m^3.s)";
      varname += "," + pp->myname + std::to_string(n + 1) + "_t";
      long_name +=
          ",horizontally averaged " + pp->myname + " " + pp->CategoryName(n);
      long_name += " density tendency";
    }
    pp = pp->next;
  }

  SetLongName(long_name);

  up_.NewAthenaArray(npart, ncells3_, ncells2_, ncells1_);
  data.NewAthenaArray(NHYDRO + npart, 1, 1, ncells1_);
  last_time_ = 0.;
}

void Tendency::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  data.ZeroClear();
  if (last_time_ == 0.) {
    wh_ = w;
    Particles *pp = pmb->ppart;
    int m = 0;
    while (pp != nullptr) {
      for (int n = 0; n < pp->u.GetDim4(); ++n)
        for (int k = ks; k <= ke; ++k)
          for (int j = js; j <= je; ++j)
            for (int i = is; i <= ie; ++i)
              up_(m + n, k, j, i) = pp->u(n, k, j, i);
      m += pp->u.GetDim4();
      pp = pp->next;
    }
    last_time_ = pmb->pmy_mesh->time;
    return;
  }

  // hydro variables
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k, j, is, ie, vol_);
      for (int n = 0; n < NHYDRO; ++n)
        for (int i = is; i <= ie; ++i)
          data(n, i) += (w(n, k, j, i) - wh_(n, k, j, i)) * vol_(i);
    }

  // particles
  Particles *pp = pmb->ppart;
  int m = 0;
  while (pp != nullptr) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        pmb->pcoord->CellVolume(k, j, is, ie, vol_);
        for (int n = 0; n < pp->u.GetDim4(); ++n)
          for (int i = is; i <= ie; ++i)
            data(NHYDRO + m + n, i) +=
                (pp->u(n, k, j, i) - up_(m + n, k, j, i)) * vol_(i);
      }
    m += pp->u.GetDim4();
    pp = pp->next;
  }

  gatherAllData23_(total_vol_, data);

  // std::cout << total_vol_(pmb->is) << std::endl;

  for (int n = 0; n < NHYDRO + m; ++n)
    for (int i = is; i <= ie; ++i)
      data(n, i) /= total_vol_(i) * (pmb->pmy_mesh->time - last_time_);

  last_time_ = pmb->pmy_mesh->time;
}
