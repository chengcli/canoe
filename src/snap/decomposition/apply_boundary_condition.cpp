// Athena++ headers
#include <coordinates/coordinates.hpp>
#include <mesh/mesh.hpp>

// canoe headers
#include "decomposition.hpp"

void Decomposition::ApplyHydroBoundary(AthenaArray<Real> &w,
                                       AthenaArray<Real> &psf, int kl, int ku,
                                       int jl, int ju) {
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Real grav = -pmy_hydro->hsrc.GetG1();  // positive downward pointing
  int is = pmb->is, ie = pmb->ie;

  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          psf(k, j, is - i) = psf(k, j, is + i);
          w(IPR, k, j, is - i) = w(IPR, k, j, is + i - 1);
          w(IDN, k, j, is - i) = -w(IDN, k, j, is + i - 1);
        }
  } else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          psf(k, j, is - i) = psf(k, j, is - i + 1) +
                              grav * w(IDN, k, j, is - i) * pco->dx1f(is - i);
          w(IPR, k, j, is - i) = w(IPR, k, j, is);
          w(IDN, k, j, is - i) = 0.;
        }
  }

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          psf(k, j, ie + i + 1) = psf(k, j, ie - i);
          w(IPR, k, j, ie + i) = w(IPR, k, j, ie - i + 1);
          w(IDN, k, j, ie + i) = -w(IDN, k, j, ie - i + 1);
        }
  } else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          psf(k, j, ie + i + 1) = psf(k, j, ie + i) - grav *
                                                          w(IDN, k, j, ie + i) *
                                                          pco->dx1f(ie + i);
          w(IPR, k, j, ie + i) = w(IPR, k, j, ie);
          w(IDN, k, j, ie + i) = 0.;
        }
  }
}

/*void Decomposition::UpdateBoundaryCondition(AthenaArray<Real> &w,
  AthenaArray<Real> const& psf, int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Real grav = -pmy_hydro->hsrc.GetG1();  // positive downward pointing

  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          w(IPR,k,j,is-i) = w(IPR,k,j,is+i-1);
          w(IDN,k,j,is-i) = w(IDN,k,j,is+i-1);
        }
  } else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          Real psv = (psf_(k,j,i) -
psf_(k,j,i+1))/log(psf_(k,j,i)/psf_(k,j,i+1)); w(IPR,k,j,is-i) += psv;
          w(IDN,k,j,is-i) = (psf_(k,j,is-i) -
psf_(k,j,is-i+1))/(pco->dx1f(is-i)*grav);
        }
  }

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          w(IPR,k,j,ie+i) = w(IPR,k,j,ie-i+1);
          w(IDN,k,j,ie+i) = w(IDN,k,j,ie-i+1);
        }
  } else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          Real psv = (psf_(k,j,ie+i) -
psf_(k,j,ie+i+1))/log(psf_(k,j,ie+i)/psf_(k,j,ie+i+1)); w(IPR,k,j,ie+i) += psv;
          w(IDN,k,j,ie+i) = (psf_(k,j,ie+i) -
psf_(k,j,ie+i+1))/(pco->dx1f(ie+i)*grav);
        }
  }
}*/
