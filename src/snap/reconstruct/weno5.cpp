// canoe
#include <configure.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
// #include "interp_weno5.hpp"
// #include "interp_weno3.hpp"

// snap
#include "face_reconstruct.hpp"
#include "interpolation.hpp"

void FaceReconstruct::Weno5X1(const int k, const int j, const int il,
                              const int iu, const AthenaArray<Real> &w,
                              const AthenaArray<Real> &bcc,
                              AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  MeshBlock *pmb = pmy_block_;
  for (int n = 0; n <= NVAPOR; ++n) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n, k, j, i + m) + 1.E-16) / 5.;
      wl(n, i + 1) =
          interp_weno5(w(n, k, j, i + 2) / scale, w(n, k, j, i + 1) / scale,
                       w(n, k, j, i) / scale, w(n, k, j, i - 1) / scale,
                       w(n, k, j, i - 2) / scale);
      wl(n, i + 1) *= scale;
      wr(n, i) =
          interp_weno5(w(n, k, j, i - 2) / scale, w(n, k, j, i - 1) / scale,
                       w(n, k, j, i) / scale, w(n, k, j, i + 1) / scale,
                       w(n, k, j, i + 2) / scale);
      wr(n, i) *= scale;
    }
  }

  int ng1 = 0, ng2 = 0;
  if (pmb->pbval->block_bcs[inner_x1] != BoundaryFlag::block) ng1 = NGHOST;
  if (pmb->pbval->block_bcs[outer_x1] != BoundaryFlag::block) ng2 = NGHOST;

  for (int n = IVX; n < NHYDRO; ++n) {
    // left boundary
    for (int i = il; i < il + ng1; ++i) {
      wl(n, i + 1) =
          interp_weno5(w(n, k, j, i + 2), w(n, k, j, i + 1), w(n, k, j, i),
                       w(n, k, j, i - 1), w(n, k, j, i - 2));
      wr(n, i) =
          interp_weno5(w(n, k, j, i - 2), w(n, k, j, i - 1), w(n, k, j, i),
                       w(n, k, j, i + 1), w(n, k, j, i + 2));
    }

    // interior
#pragma omp simd
    for (int i = il + ng1; i <= iu - ng2; ++i) {
      wl(n, i + 1) =
          interp_cp5(w(n, k, j, i + 2), w(n, k, j, i + 1), w(n, k, j, i),
                     w(n, k, j, i - 1), w(n, k, j, i - 2));
      wr(n, i) = interp_cp5(w(n, k, j, i - 2), w(n, k, j, i - 1), w(n, k, j, i),
                            w(n, k, j, i + 1), w(n, k, j, i + 2));
    }

    // right boundary
    for (int i = iu - ng2 + 1; i <= iu; ++i) {
      wl(n, i + 1) =
          interp_weno5(w(n, k, j, i + 2), w(n, k, j, i + 1), w(n, k, j, i),
                       w(n, k, j, i - 1), w(n, k, j, i - 2));
      wr(n, i) =
          interp_weno5(w(n, k, j, i - 2), w(n, k, j, i - 1), w(n, k, j, i),
                       w(n, k, j, i + 1), w(n, k, j, i + 2));
    }
  }

  return;
}

void FaceReconstruct::Weno5X2(const int k, const int j, const int il,
                              const int iu, const AthenaArray<Real> &w,
                              const AthenaArray<Real> &bcc,
                              AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  for (int n = 0; n <= NVAPOR; ++n) {
    for (int i = il; i <= iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n, k, j + m, i) + 1.E-16) / 5.;
      wl(n, i) =
          interp_weno5(w(n, k, j + 2, i) / scale, w(n, k, j + 1, i) / scale,
                       w(n, k, j, i) / scale, w(n, k, j - 1, i) / scale,
                       w(n, k, j - 2, i) / scale);
      wl(n, i) *= scale;
      wr(n, i) =
          interp_weno5(w(n, k, j - 2, i) / scale, w(n, k, j - 1, i) / scale,
                       w(n, k, j, i) / scale, w(n, k, j + 1, i) / scale,
                       w(n, k, j + 2, i) / scale);
      wr(n, i) *= scale;
    }
  }

  for (int n = IVX; n < NHYDRO; ++n) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      wl(n, i) = interp_cp5(w(n, k, j + 2, i), w(n, k, j + 1, i), w(n, k, j, i),
                            w(n, k, j - 1, i), w(n, k, j - 2, i));
      wr(n, i) = interp_cp5(w(n, k, j - 2, i), w(n, k, j - 1, i), w(n, k, j, i),
                            w(n, k, j + 1, i), w(n, k, j + 2, i));
    }
  }

  return;
}

void FaceReconstruct::Weno5X3(const int k, const int j, const int il,
                              const int iu, const AthenaArray<Real> &w,
                              const AthenaArray<Real> &bcc,
                              AthenaArray<Real> &wl, AthenaArray<Real> &wr) {
  for (int n = 0; n <= NVAPOR; ++n) {
    for (int i = il; i <= iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n, k + m, j, i) + 1.E-16) / 5.;
      wl(n, i) =
          interp_weno5(w(n, k + 2, j, i) / scale, w(n, k + 1, j, i) / scale,
                       w(n, k, j, i) / scale, w(n, k - 1, j, i) / scale,
                       w(n, k - 2, j, i) / scale);
      wl(n, i) *= scale;
      wr(n, i) =
          interp_weno5(w(n, k - 2, j, i) / scale, w(n, k - 1, j, i) / scale,
                       w(n, k, j, i) / scale, w(n, k + 1, j, i) / scale,
                       w(n, k + 2, j, i) / scale);
      wr(n, i) *= scale;
    }
  }

  for (int n = IVX; n < NHYDRO; ++n) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      wl(n, i) = interp_cp5(w(n, k + 2, j, i), w(n, k + 1, j, i), w(n, k, j, i),
                            w(n, k - 1, j, i), w(n, k - 2, j, i));
      wr(n, i) = interp_cp5(w(n, k - 2, j, i), w(n, k - 1, j, i), w(n, k, j, i),
                            w(n, k + 1, j, i), w(n, k + 2, j, i));
    }
  }

  return;
}
