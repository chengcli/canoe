// application
#include <application/application.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/hydro/srcterms/hydro_srcterms.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/stride_iterator.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// canoe
#include <configure.h>

#include <impl.hpp>
#include <interface/eos.hpp>

void check_hydro_variables(MeshBlock *pmb) {
  Application::Logger app("snap");

  auto &w = pmb->phydro->w;
  auto temp = get_temp(pmb->pimpl->peos, w);
  auto tempa = temp.accessor<Real, 3>();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        for (int n = 0; n <= NHYDRO; ++n) {
          if (w(n, k, j, i) < 0.) {
            app->Error("density is negative at (" + std::to_string(k) + ", " +
                       std::to_string(j) + ", " + std::to_string(i) + ")");
          }
        }

        if (NON_BAROTROPIC_EOS) {
          if (w(IPR, k, j, i) < 0.) {
            app->Error("pressure is negative at (" + std::to_string(k) + ", " +
                       std::to_string(j) + ", " + std::to_string(i) + ")");
          }
        }

        Real grav = -pmb->phydro->hsrc.GetG1();
        Real Rd = get_rd();
        if (grav != 0) {
          Real Tmin = grav * pmb->pcoord->dx1f(i) / Rd;
          if (tempa[k][j][i] < Tmin) {
            app->Error("temperature is less than minimum temperature at (" +
                       std::to_string(k) + ", " + std::to_string(j) + ", " +
                       std::to_string(i) + ")");
          }
        }
      }

  // make a copy of w, needed for outflow boundary condition
  // w1 = w;
  app->Log("Hydro check passed.");
}
