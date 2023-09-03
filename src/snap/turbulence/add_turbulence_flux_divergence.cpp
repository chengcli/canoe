// Athena++ headers
#include <athena/athena.hpp>
#include <athena/athena_arrays.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

// snap
#include "turbulence_model.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void TurbulenceModel::AddFluxDivergence
//  \brief Adds flux divergence to weighted average of conservative variables
//  from previous step(s) of time integrator algorithm

// TODO(felker): after the removal of AddCoordTermsDivergence() fn call from
// Hydro::AddFluxDivergence(), the 2x fns could be trivially shared if:
// - flux/s_flux renamed to the same class member name
// - 7x below references of x1face_area_ ... dflx_ private class members (which
// are only ever used in this fn and are members to prevent de/allocating each
// fn call)
// - NHYDRO/NSCALARS is replaced with array_out.GetDim4()

// ----> Hydro should be derived from Turbulence

// TODO(felker): remove the following unnecessary private class member?
// field_diffusion.cpp:66:    cell_volume_.NewAthenaArray(nc1);

void TurbulenceModel::addFluxDivergence(const Real wght,
                                        AthenaArray<Real> &s_out) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = s_flux[X1DIR];
  AthenaArray<Real> &x2flux = s_flux[X2DIR];
  AthenaArray<Real> &x3flux = s_flux[X3DIR];
  int is = pmb->is;
  int js = pmb->js;
  int ks = pmb->ks;
  int ie = pmb->ie;
  int je = pmb->je;
  int ke = pmb->ke;
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                    &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                    &x3area_p1 = x3face_area_p1_, &vol = cell_volume_,
                    &dflx = dflx_;

  int nvar = s_out.GetDim4();
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, is, ie + 1, x1area);
      for (int n = 0; n < nvar; ++n) {
#pragma omp simd
        for (int i = is; i <= ie; ++i) {
          dflx(n, i) = (x1area(i + 1) * x1flux(n, k, j, i + 1) -
                        x1area(i) * x1flux(n, k, j, i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j, is, ie, x2area);
        pmb->pcoord->Face2Area(k, j + 1, is, ie, x2area_p1);
        for (int n = 0; n < nvar; ++n) {
#pragma omp simd
          for (int i = is; i <= ie; ++i) {
            dflx(n, i) += (x2area_p1(i) * x2flux(n, k, j + 1, i) -
                           x2area(i) * x2flux(n, k, j, i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k, j, is, ie, x3area);
        pmb->pcoord->Face3Area(k + 1, j, is, ie, x3area_p1);
        for (int n = 0; n < nvar; ++n) {
#pragma omp simd
          for (int i = is; i <= ie; ++i) {
            dflx(n, i) += (x3area_p1(i) * x3flux(n, k + 1, j, i) -
                           x3area(i) * x3flux(n, k, j, i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k, j, is, ie, vol);
      for (int n = 0; n < nvar; ++n) {
#pragma omp simd
        for (int i = is; i <= ie; ++i) {
          s_out(n, k, j, i) -= wght * dflx(n, i) / vol(i);
        }
      }
    }
  }
  return;
}
