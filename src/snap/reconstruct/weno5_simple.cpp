//! \file weno5_simple.cpp
//  \brief  WENO5 interpolation

// Athena++ headers
#include <athena.hpp>
#include <mesh/mesh.hpp>

// canoe headers
#include <configure.hpp>

#include "face_reconstruct.hpp"
#include "interpolation.hpp"

void FaceReconstruct::Weno5X1Simple(const int k, const int j, const int il,
                                    const int iu, const AthenaArray<Real> &w,
                                    AthenaArray<Real> &wl,
                                    AthenaArray<Real> &wr) {
  for (int n = 0; n < w.GetDim4(); ++n)
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      wl(n, i + 1) =
          interp_weno5(w(n, k, j, i + 2), w(n, k, j, i + 1), w(n, k, j, i),
                       w(n, k, j, i - 1), w(n, k, j, i - 2));
      wr(n, i) =
          interp_weno5(w(n, k, j, i - 2), w(n, k, j, i - 1), w(n, k, j, i),
                       w(n, k, j, i + 1), w(n, k, j, i + 2));
    }

  return;
}

void FaceReconstruct::Weno5X2Simple(const int k, const int j, const int il,
                                    const int iu, const AthenaArray<Real> &w,
                                    AthenaArray<Real> &wl,
                                    AthenaArray<Real> &wr) {
  for (int n = 0; n < w.GetDim4(); ++n)
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      wl(n, i) =
          interp_weno5(w(n, k, j + 2, i), w(n, k, j + 1, i), w(n, k, j, i),
                       w(n, k, j - 1, i), w(n, k, j - 2, i));
      wr(n, i) =
          interp_weno5(w(n, k, j - 2, i), w(n, k, j - 1, i), w(n, k, j, i),
                       w(n, k, j + 1, i), w(n, k, j + 2, i));
    }

  return;
}

void FaceReconstruct::Weno5X3Simple(const int k, const int j, const int il,
                                    const int iu, const AthenaArray<Real> &w,
                                    AthenaArray<Real> &wl,
                                    AthenaArray<Real> &wr) {
  for (int n = 0; n < w.GetDim4(); ++n)
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      wl(n, i) =
          interp_weno5(w(n, k + 2, j, i), w(n, k + 1, j, i), w(n, k, j, i),
                       w(n, k - 1, j, i), w(n, k - 2, j, i));
      wr(n, i) =
          interp_weno5(w(n, k - 2, j, i), w(n, k - 1, j, i), w(n, k, j, i),
                       w(n, k + 1, j, i), w(n, k + 2, j, i));
    }

  return;
}
