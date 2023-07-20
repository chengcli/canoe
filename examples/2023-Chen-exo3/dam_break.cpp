// C/C++
#include <cmath>

// Athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <configure.hpp>

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real phi = pin->GetReal("problem", "phi");
  // Real uphi = pin->GetReal("problem", "uphi");
  Real dphi = pin->GetReal("problem", "dphi");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      phydro->u(IDN, k, j, is) = phi;
      // if (pcoord->x3v(k) > 0.)
      //   phydro->u(IM2,k,j,is) = -uphi;
      // else
      //   phydro->u(IM2,k,j,is) = uphi;
#ifdef AFFINE
      Real x1_car = pcoord->x2v(j) + pcoord->x3v(k) * cos(M_PI / 3.0);
      Real x2_car = pcoord->x3v(k) * sin(M_PI / 3.0);
      Real rad = sqrt(x1_car * x1_car + x2_car * x2_car);
#else
      Real rad = sqrt(pcoord->x3v(k) * pcoord->x3v(k) +
                      pcoord->x2v(j) * pcoord->x2v(j));
#endif
      if (rad < 2.)  // Circular dam breaking
        phydro->u(IDN, k, j, is) += dphi;
    }
}
