/* Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "interpolation.hpp"
#include "reconstruction.hpp"

using namespace torch;

// C,D,H,W
enum {
  DIMC = 0,
  DIM3 = 1,
  DIM2 = 2,
  DIM1 = 3,
}

void Reconstruction::Weno5X1(int il, int iu,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  MeshBlock *pmb = pmy_block_;

  Tensor *w_ = w.Tensor(DIMC, {0, IVX});
  Tensor *wl_ = wl.Tensor(DIMC, {0, IVX});
  Tensor *wr_ = wr.Tensor(DIMC, {0, IVX});

  for (int i=il; i<=iu; ++i) {
    auto &wi = w_->slice(DIM1, i-2, i+2);
    Tensor scale = wi.sum(DIM1) + 1.E-6;
    wi /= scale.unsqueeze(DIM1);

    wl_->narrow(DIM1, i+1, 1) = interp_weno5p(wi, "nkji,i->nkj");
    wl_->narrow(DIM1, i+1, 1) *= scale;

    wr_->narrow(DIM1, i, 1) = interp_weno5m(wi, "nkji,i->nkj");
    wr_->narrow(DIM1, i, 1) *= scale;

    wi *= scale.unsqueeze(DIM1);
  };

  if (shock_capture_flag_) {
    w_ = w.Tensor(DIMC, {IVX, NHYDRO});
    wl_ = wl.Tensor(DIMC, {IVX, NHYDRO});
    wr_ = wr.Tensor(DIMC, {IVX, NHYDRO});
  }

    wl(n,i+1) = interp_weno5(w(n,k,j,i+2)/scale, w(n,k,j,i+1)/scale,
w(n,k,j,i)/scale, w(n,k,j,i-1)/scale, w(n,k,j,i-2)/scale); wl(n,i+1) *= scale;
    wr(n,i) = interp_weno5(w(n,k,j,i-2)/scale, w(n,k,j,i-1)/scale,
w(n,k,j,i)/scale, w(n,k,j,i+1)/scale, w(n,k,j,i+2)/scale); wr(n,i) *= scale;
  }

  if (shock_capture_flag_) {
    for (int n=IVX; n<NHYDRO; ++n)
      for (int i=il; i<=iu; ++i) {
        wl(n,i+1) =
interp_weno5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
        wr(n,i  ) =
interp_weno5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
      }
    return;
  }

  int ng1 = 0, ng2 = 0;
  if (pmb->pbval->isPhysicalBoundary(inner_x1))
    ng1 = NGHOST;
  if (pmb->pbval->isPhysicalBoundary(outer_x1))
    ng2 = NGHOST;

  for (int n=IVX; n<NHYDRO; ++n) {
    // left boundary
    for (int i=il; i<il+ng1; ++i) {
      wl(n,i+1) =
interp_weno5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i  ) =
interp_weno5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }

    // interior
#pragma omp simd
    for (int i=il+ng1; i<=iu-ng2; ++i) {
      wl(n,i+1) =
interp_cp5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i  ) =
interp_cp5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }

    // right boundary
    for (int i=iu-ng2+1; i<=iu; ++i) {
      wl(n,i+1) =
interp_weno5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i  ) =
interp_weno5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X2()
//  \brief

void Reconstruction::Weno5X2(int jl, int ju,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  auto w_ = w.Tensor(DIMC, {0, IVX});
  auto wl_ = wl.Tensor(DIMC, {0, IVX});
  auto wr_ = wr.Tensor(DIMC, {0, IVX});

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

  for (int n=0; n<IVX; ++n) {
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n,k,j+m,i) + 1.E-16)/5.;
      wl(n,i) = interp_weno5(w(n,k,j+2,i)/scale, w(n,k,j+1,i)/scale,
w(n,k,j,i)/scale, w(n,k,j-1,i)/scale, w(n,k,j-2,i)/scale); wl(n,i) *= scale;
      wr(n,i) = interp_weno5(w(n,k,j-2,i)/scale, w(n,k,j-1,i)/scale,
w(n,k,j,i)/scale, w(n,k,j+1,i)/scale, w(n,k,j+2,i)/scale); wr(n,i) *= scale;
    }
  }

  if (shock_capture_flag_) {
    for (int n=IVX; n<NHYDRO; ++n)
      for (int i=il; i<=iu; ++i) {
        wl(n,i) =
interp_weno5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
        wr(n,i) =
interp_weno5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
      }
    return;
  }

  for (int n=IVX; n<NHYDRO; ++n) {
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
