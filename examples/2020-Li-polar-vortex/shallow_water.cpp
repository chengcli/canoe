//! \file shallow_water.cpp
//  \brief shallow water test model

// C/C++
#include <cmath>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real phi = pin->GetReal("problem", "phi");
  Real uphi = pin->GetReal("problem", "uphi");
  Real dphi = pin->GetReal("problem", "dphi");

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN, k, j, i) = phi;
        if (pcoord->x2v(j) > 0.)
          phydro->u(IM1, k, j, i) = -uphi;
        else
          phydro->u(IM1, k, j, i) = uphi;

        if (pcoord->x1v(i) > 0. && pcoord->x1v(i) < 5.)
          phydro->u(IDN, k, j, i) += dphi;
      }
}
