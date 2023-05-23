// Athena++ headers
#include <athena.hpp>
#include <coordinates/coordinates.hpp>
#include <hydro/hydro.hpp>
#include <hydro/srcterms/hydro_srcterms.hpp>
#include <mesh/mesh.hpp>

// debugger headers
#include <debugger.hpp>

// cliutils headers
#include <cliutils/stride_iterator.hpp>

// canoe headers
#include <configure.hpp>

#include "../mesh/meshblock_impl.hpp"
#include "../thermodynamics/thermodynamics.hpp"

void check_hydro_variables(MeshBlock *pmb, AthenaArray<Real> const &w) {
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        for (int n = 0; n <= NumVapors; ++n) {
          if (w(n, k, j, i) < 0.) {
            Debugger::Fatal("check_hydro_variables", "density", "negative");
          }
        }

#if NON_BAROTROPIC_EOS == 1
        if (w(IPR, k, j, i) < 0.) {
          Debugger::Fatal("check_hydro_variables", "pressure", "negative");
        }
#endif

        Thermodynamics *pthermo = pmb->pimpl->pthermo;
        Real temp = pthermo->getTemp(w.at(k, j, i));
        Real grav = -pmb->phydro->hsrc.GetG1();
        if (grav != 0) {
          Real Tmin = grav * pmb->pcoord->dx1f(i) / pthermo->getRd();
          if (temp < Tmin) {
            Debugger::Fatal("check_hydro_variables", "temperature", "low");
          }
        }
      }

  // make a copy of w, needed for outflow boundary condition
  // w1 = w;
  Debugger::Print("Hydro check passed.");
}
