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
/*  \brief

void Reconstruction::Weno5X2(int jl, int ju,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  int64_t nweno = shock_capture_flag_ ? IVX : NHYDRO;

  auto w_ = w.Tensor(DIMC, {0, nweno});
  auto wl_ = wl.Tensor(DIMC, {0, nweno});
  auto wr_ = wr.Tensor(DIMC, {0, nweno});

  for (int j=jl; j<=ju; ++j) {
    auto &wj = w_->slice(DIM2, j-2, j+2);
    Tensor scale = wj.sum(DIM2) + 1.E-6;

    wj /= scale.unsqueeze(DIM2);
    wl_->narrow(DIM2, j, 1) = interp_weno5(wj, "nkji,j->nki");
    wl_->narrow(DIM2, j, 1) *= scale;

    wr_->narrow(DIM2, j, 1) = interp_weno5p(wj, "nkji,j->nki");
    wr_->narrow(DIM2, j, 1) *= scale;

    wj *= scale.unsqueeze(DIM2);
  }

  for (int n=nweno; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) =
interp_cp5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
      wr(n,i) =
interp_cp5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X3()
//  \brief

void Reconstruction::Weno5X3(int kl, int ku,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  auto w_ = w.Tensor(DIMC, {0, IVX});
  auto wl_ = wl.Tensor(DIMC, {0, IVX});
  auto wr_ = wr.Tensor(DIMC, {0, IVX});

  for (int k=kl; k<=ku; ++k) {
    auto &wk = w_->slice(DIM3, k-2, k+2);
    Tensor scale = wk.sum() + 1.E-6;
    wk /= scale.unsqueeze(DIM3);

    wl_->narrow(DIM3, k, 1) = interp_weno5(wk, "nkji,k->nji");
    wr_->narrow(DIM3, k, 1) *= scale;

    wr_->narrow(DIM3, k, 1) = interp_weno5p(wk, "nkji,k->nji");
    wr_->narrow(DIM3, k, 1) *= scale;

    wk *= scale.unsqueeze(DIM3);
  }

  for (int n=0; n<IVX; ++n) {
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n,k+m,j,i) + 1.E-16)/5.;
      wl(n,i) = interp_weno5(w(n,k+2,j,i)/scale, w(n,k+1,j,i)/scale,
w(n,k,j,i)/scale, w(n,k-1,j,i)/scale, w(n,k-2,j,i)/scale); wl(n,i) *= scale;
      wr(n,i) = interp_weno5(w(n,k-2,j,i)/scale, w(n,k-1,j,i)/scale,
w(n,k,j,i)/scale, w(n,k+1,j,i)/scale, w(n,k+2,j,i)/scale); wr(n,i) *= scale;
    }
  }

  if (shock_capture_flag_) {
    for (int n=IVX; n<NHYDRO; ++n)
      for (int i=il; i<=iu; ++i) {
        wl(n,i) =
interp_weno5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
        wr(n,i) =
interp_weno5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
      }
    return;
  }

  for (int n=IVX; n<NHYDRO; ++n) {
    for (int i=il; i<=iu; ++i) {
      wl(n,i) =
interp_cp5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
      wr(n,i) =
interp_cp5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
    }
  }
}*/
