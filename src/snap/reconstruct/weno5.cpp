// C/C++ header
#include <limits>

// snap
#include "interpolation.hpp"

using namespace torch;

// C,D,H,W
enum {
  DIMC = 0,
  DIM3 = 1,
  DIM2 = 2,
  DIM1 = 3,
  DIMX = 4,
  NGHOST = 3,
};

inline void _weno5_hydro_range(int64_t dim, int64_t il, int64_t iu,
                               const Tensor &w, Tensor &wl, Tensor &wr,
                               const Interpolation &interp) {
  if (il > iu) return;

  auto wu = w.slice(dim, il - NGHOST, iu + NGHOST - 1).unfold(dim, 5, 1);
  auto scale = wu.abs().mean(DIMX) + std::numeric_limits<float>::min();
  wu /= scale.unsqueeze(DIMX);

  wl.slice(dim, il, iu + 1) = interp.right(wu) * scale;
  wr.slice(dim, il - 1, iu) = interp.left(wu) * scale;
}

std::pair<Tensor, Tensor> recon_weno5_hydro(const Tensor &w, int64_t IVX,
                                            int64_t dim, bool mixed,
                                            bool is_boundary_lower,
                                            bool is_boundary_upper) {
  int64_t NHYDRO = w.size(DIMC);
  int64_t nweno = mixed ? IVX : NHYDRO;
  Weno5Interp interp1(w.device().type());
  Center5Interp interp2(w.device().type());

  auto wl = torch::zeros_like(w);
  auto wr = torch::zeros_like(w);

  auto w_ = w.slice(DIMC, 0, nweno);
  auto wl_ = wl.slice(DIMC, 0, nweno);
  auto wr_ = wr.slice(DIMC, 0, nweno);

  auto dim_size = w.size(dim);
  int64_t il = NGHOST;
  int64_t iu = dim_size - NGHOST + 1;

  _weno5_hydro_range(dim, il, iu, w_, wl_, wr_, interp1);

  if (nweno == NHYDRO) return {wl, wr};

  // rest of the hydro variables
  w_ = w.slice(DIMC, nweno, NHYDRO);
  wl_ = wl.slice(DIMC, nweno, NHYDRO);
  wr_ = wr.slice(DIMC, nweno, NHYDRO);

  // boundaries
  if (dim_size > 2 * NGHOST) {
    if (is_boundary_lower) {
      _weno5_hydro_range(dim, il, il + NGHOST - 1, w_, wl_, wr_, interp1);
      il += NGHOST;
    } else if (is_boundary_upper) {
      _weno5_hydro_range(dim, iu - NGHOST + 1, iu, w_, wl_, wr_, interp1);
      iu -= NGHOST;
    }
  } else {
    if (is_boundary_lower && !is_boundary_upper) {
      _weno5_hydro_range(dim, il, il + NGHOST - 1, w_, wl_, wr_, interp1);
      il += NGHOST;
    } else if (!is_boundary_lower && is_boundary_upper) {
      _weno5_hydro_range(dim, iu - NGHOST + 1, iu, w_, wl_, wr_, interp1);
      iu -= NGHOST;
    } else if (is_boundary_lower && is_boundary_upper) {
      int64_t len1 = dim_size / 2;
      int64_t len2 = dim_size - len1;
      _weno5_hydro_range(dim, il, il + len1 - 1, w_, wl_, wr_, interp1);
      _weno5_hydro_range(dim, iu - len2 + 1, iu, w_, wl_, wr_, interp1);
      il += len1;
      iu -= len2;
    }
  }

  // interior
  _weno5_hydro_range(dim, il, iu, w_, wl_, wr_, interp2);

  return {wl, wr};
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X2()
/*  \brief

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
}*/
