// C/C++ header
#include <limits>

// Athena++ headers
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/reconstruct/reconstruction.hpp>

// snap
#include <snap/athena_arrays.hpp>

#include "interpolation.hpp"

using namespace torch;

// C,D,H,W
enum {
  DIMC = 0,
  DIM3 = 1,
  DIM2 = 2,
  DIM1 = 3,
};

void Reconstruction::Weno5X1(int il, int iu, const Tensor &w, Tensor &wl,
                             Tensor &wr) {
  int64_t nweno = shock_capture_flag_ ? NHYDRO : IVX;
  Weno5Interp interp1(w.device().type());
  Center5Interp interp2(w.device().type());

  auto w_ = w.slice(DIMC, 0, nweno);
  auto wl_ = wl.slice(DIMC, 0, nweno);
  auto wr_ = wr.slice(DIMC, 0, nweno);

  for (int i = il; i <= iu; ++i) {
    auto wi = w_.slice(DIM1, i - 2, i + 3);
    auto scale = wi.abs().mean(DIM1) + std::numeric_limits<float>::min();
    scale = scale.unsqueeze(DIM1);
    wi /= scale;

    wl_.narrow(DIM1, i + 1, 1) =
        interp1.right(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wl_.narrow(DIM1, i + 1, 1) *= scale;

    wr_.narrow(DIM1, i, 1) = interp1.left(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wr_.narrow(DIM1, i, 1) *= scale;

    wi *= scale;
  };

  MeshBlock *pmb = pmy_block_;
  int ng1 = pmb->pbval->isPhysicalBoundary(inner_x1) ? NGHOST : 0;
  int ng2 = pmb->pbval->isPhysicalBoundary(outer_x1) ? NGHOST : 0;

  if (nweno == NHYDRO) return;

  // rest of the hydro variables
  w_ = w.slice(DIMC, nweno, NHYDRO);
  wl_ = wl.slice(DIMC, nweno, NHYDRO);
  wr_ = wr.slice(DIMC, nweno, NHYDRO);

  // left boundary
  for (int i = il; i < il + ng1; ++i) {
    auto wi = w_.slice(DIM1, i - 2, i + 3);
    auto scale = wi.abs().mean(DIM1) + std::numeric_limits<float>::min();
    scale = scale.unsqueeze(DIM1);
    wi /= scale;
    wl_.narrow(DIM1, i + 1, 1) =
        interp1.right(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wl_.narrow(DIM1, i + 1, 1) *= scale;
    wr_.narrow(DIM1, i, 1) = interp1.left(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wr_.narrow(DIM1, i, 1) *= scale;
    wi *= scale;
  }

  // interior
  for (int i = il + ng1; i <= iu - ng2; ++i) {
    auto wi = w_.slice(DIM1, i - 2, i + 3);
    wl_.narrow(DIM1, i + 1, 1) =
        interp2.right(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wr_.narrow(DIM1, i, 1) = interp2.left(wi, "nkji,i->nkj").unsqueeze(DIM1);
  }

  // right boundary
  for (int i = iu - ng2 + 1; i <= iu; ++i) {
    auto wi = w_.slice(DIM1, i - 2, i + 3);
    auto scale = wi.abs().mean(DIM1) + std::numeric_limits<float>::min();
    scale = scale.unsqueeze(DIM1);
    wi /= scale;
    wl_.narrow(DIM1, i + 1, 1) =
        interp1.right(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wl_.narrow(DIM1, i + 1, 1) *= scale;
    wr_.narrow(DIM1, i, 1) = interp1.left(wi, "nkji,i->nkj").unsqueeze(DIM1);
    wr_.narrow(DIM1, i, 1) *= scale;
    wi *= scale;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X2()
//  \brief

void Reconstruction::Weno5X2(int jl, int ju, const Tensor &w, Tensor &wl,
                             Tensor &wr) {
  int64_t nweno = shock_capture_flag_ ? NHYDRO : IVX;
  Weno5Interp interp1(w.device().type());
  Center5Interp interp2(w.device().type());

  auto w_ = w.slice(DIMC, 0, nweno);
  auto wl_ = wl.slice(DIMC, 0, nweno);
  auto wr_ = wr.slice(DIMC, 0, nweno);

  for (int j = jl; j <= ju; ++j) {
    auto wj = w_.slice(DIM2, j - 2, j + 3);
    auto scale = wj.abs().mean(DIM2) + std::numeric_limits<float>::min();
    scale = scale.unsqueeze(DIM2);
    wj /= scale;

    wl_.narrow(DIM2, j + 1, 1) =
        interp1.right(wj, "nkji,j->nki").unsqueeze(DIM2);
    wl_.narrow(DIM2, j + 1, 1) *= scale;

    wr_.narrow(DIM2, j, 1) = interp1.left(wj, "nkji,j->nki").unsqueeze(DIM2);
    wr_.narrow(DIM2, j, 1) *= scale;

    wj *= scale;
  }

  if (nweno == NHYDRO) return;

  // rest of the hydro variables
  w_ = w.slice(DIMC, nweno, NHYDRO);
  wl_ = wl.slice(DIMC, nweno, NHYDRO);
  wr_ = wr.slice(DIMC, nweno, NHYDRO);

  for (int j = jl; j <= ju; ++j) {
    auto wj = w_.slice(DIM2, j - 2, j + 3);
    wl_.narrow(DIM2, j + 1, 1) =
        interp2.right(wj, "nkji,j->nki").unsqueeze(DIM2);
    wr_.narrow(DIM2, j, 1) = interp2.left(wj, "nkji,j->nki").unsqueeze(DIM2);
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X3()
//  \brief

void Reconstruction::Weno5X3(int kl, int ku, const Tensor &w, Tensor &wl,
                             Tensor &wr) {
  int64_t nweno = shock_capture_flag_ ? NHYDRO : IVX;
  Weno5Interp interp1(w.device().type());
  Center5Interp interp2(w.device().type());

  auto w_ = w.slice(DIMC, 0, nweno);
  auto wl_ = wl.slice(DIMC, 0, nweno);
  auto wr_ = wr.slice(DIMC, 0, nweno);

  for (int k = kl; k <= ku; ++k) {
    auto wk = w_.slice(DIM3, k - 2, k + 3);
    auto scale = wk.abs().mean(DIM3) + std::numeric_limits<float>::min();
    scale = scale.unsqueeze(DIM3);
    wk /= scale;

    wl_.narrow(DIM3, k + 1, 1) =
        interp1.right(wk, "nkji,k->nji").unsqueeze(DIM3);
    wr_.narrow(DIM3, k + 1, 1) *= scale;

    wr_.narrow(DIM3, k, 1) = interp1.left(wk, "nkji,k->nji").unsqueeze(DIM3);
    wr_.narrow(DIM3, k, 1) *= scale;

    wk *= scale;
  }

  if (nweno == NHYDRO) return;

  // rest of the hydro variables
  w_ = w.slice(DIMC, nweno, NHYDRO);
  wl_ = wl.slice(DIMC, nweno, NHYDRO);
  wr_ = wr.slice(DIMC, nweno, NHYDRO);

  for (int k = kl; k <= ku; ++k) {
    auto wk = w_.slice(DIM3, k - 2, k + 3);
    wl_.narrow(DIM3, k + 1, 1) =
        interp2.right(wk, "nkji,k->nji").unsqueeze(DIM3);
    wr_.narrow(DIM3, k, 1) = interp2.left(wk, "nkji,k->nji").unsqueeze(DIM3);
  }
}
