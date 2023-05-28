// dealii headers
#include <deal.II/lac/precondition.h>

// athena headers
#include <athena/athena.hpp>

// IceChem headers
#include "stream_transport.hpp"

StreamTransport::StreamTransport(int rows, int cols, bool fourth_order)
    : fourth_order_(fourth_order),
      rows_(rows),
      cols_(cols),
      rank_(rows * cols),
      rowsh_(rows + 2 * GhostZoneSize),
      colsh_(cols + 2 * GhostZoneSize),
      rankh_(rowsh_ * colsh_),
      skh_(rank_, rankh_, 9),
      skk_(rank_, rank_, 9),
      shk_(rankh_, rank_, 1),
      qvec_(rankh_),
      dqvec_(rank_),
      rhs_(rank_),
      control_(1000, 1e-12),
      solver_(control_) {
  // 9 elements maximum
  for (int j = 0; j < rows_; ++j)
    for (int i = 0; i < cols_; ++i) {
      int k = global(j, i);
      skh_.add(k, globalh(j, i));
      skh_.add(k, globalh(j - 1, i));
      skh_.add(k, globalh(j, i - 1));
      skh_.add(k, globalh(j + 1, i));
      skh_.add(k, globalh(j, i + 1));
      skh_.add(k, globalh(j - 1, i - 1));
      skh_.add(k, globalh(j - 1, i + 1));
      skh_.add(k, globalh(j + 1, i - 1));
      skh_.add(k, globalh(j + 1, i + 1));
    }
  skh_.compress();

  diffusion_.reinit(skh_);
  advection_.reinit(skh_);
  KmJ_.reinit(skh_);

  // 9 elements maximum
  for (int j = 0; j < rows_; ++j)
    for (int i = 0; i < cols_; ++i) {
      int64_t k = global(j, i);
      skk_.add(k, global(j, i));
      if (j > 0) skk_.add(k, global(j - 1, i));
      if (i > 0) skk_.add(k, global(j, i - 1));
      if (j < rows_ - 1) skk_.add(k, global(j + 1, i));
      if (i < cols_ - 1) skk_.add(k, global(j, i + 1));
      if (j > 0 && i > 0) skk_.add(k, global(j - 1, i - 1));
      if (j > 0 && i < cols_ - 1) skk_.add(k, global(j - 1, i + 1));
      if (j < rows_ - 1 && i > 0) skk_.add(k, global(j + 1, i - 1));
      if (j < rows_ - 1 && i < cols_ - 1) skk_.add(k, global(j + 1, i + 1));
    }
  skk_.compress();

  mass_.reinit(skk_);
  KmJmN_.reinit(skk_);

  // 2 element
  for (int j = 0; j < rows_; ++j)
    for (int i = 0; i < cols_; ++i) {
      shk_.add(globalh(j, i), global(j, i));
    }
  shk_.compress();

  // neumann boundary condition
  bneumann_.reinit(shk_);
  for (int i = 0; i < cols_; ++i) {
    bneumann_.set(globalh(-1, i), global(0, i), 1.);
    for (int j = 0; j < rows_; ++j)
      bneumann_.set(globalh(j, i), global(j, i), 1.);
    bneumann_.set(globalh(rows_, i), global(rows_ - 1, i), 1.);
  }
}

StreamTransport::~StreamTransport() {}

void StreamTransport::setDiffusionMatrix(Real kdiff, Real dx) {
  if (fourth_order_) {
    for (int j = 0; j < rows_; ++j)
      for (int i = 0; i < cols_; ++i) {
        int64_t k = global(j, i);
        diffusion_.set(k, globalh(j, i), -10. / 3.);
        diffusion_.set(k, globalh(j - 1, i), 2. / 3.);
        diffusion_.set(k, globalh(j, i - 1), 2. / 3.);
        diffusion_.set(k, globalh(j + 1, i), 2. / 3.);
        diffusion_.set(k, globalh(j, i + 1), 2. / 3.);
        diffusion_.set(k, globalh(j - 1, i - 1), 1. / 6.);
        diffusion_.set(k, globalh(j - 1, i + 1), 1. / 6.);
        diffusion_.set(k, globalh(j + 1, i - 1), 1. / 6.);
        diffusion_.set(k, globalh(j + 1, i + 1), 1. / 6.);
      }
  } else {
    for (int j = 0; j < rows_; ++j)
      for (int i = 0; i < cols_; ++i) {
        int64_t k = global(j, i);
        diffusion_.set(k, globalh(j, i), -4.);
        diffusion_.set(k, globalh(j - 1, i), 1.);
        diffusion_.set(k, globalh(j + 1, i), 1.);
        diffusion_.set(k, globalh(j, i - 1), 1.);
        diffusion_.set(k, globalh(j, i + 1), 1.);
      }
  }
  diffusion_ *= kdiff / (dx * dx);
}

void StreamTransport::setAdvectionMatrix(AthenaArray<Real> const& streamf,
                                         Real dx) {
  for (int j = 0; j < rows_; ++j)
    for (int i = 0; i < cols_; ++i) {
      int64_t k = global(j, i);
      int j1 = j + GhostZoneSize;
      int i1 = i + GhostZoneSize;
      advection_.set(k, globalh(j - 1, i - 1),
                     -streamf(j1, i1 - 1) + streamf(j1 - 1, i1));

      advection_.set(k, globalh(j - 1, i),
                     -streamf(j1 - 1, i1 - 1) - streamf(j1, i1 - 1) +
                         streamf(j1 - 1, i1 + 1) + streamf(j1, i1 + 1));

      advection_.set(k, globalh(j - 1, i + 1),
                     +streamf(j1, i1 + 1) - streamf(j1 - 1, i1));

      advection_.set(k, globalh(j, i - 1),
                     -streamf(j1 + 1, i1 - 1) - streamf(j1 + 1, i1) +
                         streamf(j1 - 1, i1 - 1) + streamf(j1 - 1, i1));

      advection_.set(k, globalh(j, i + 1),
                     +streamf(j1 + 1, i1) + streamf(j1 + 1, i1 + 1) -
                         streamf(j1 - 1, i1) - streamf(j1 - 1, i1 + 1));

      advection_.set(k, globalh(j + 1, i - 1),
                     -streamf(j1 + 1, i1) + streamf(j1, i1 - 1));

      advection_.set(k, globalh(j + 1, i),
                     +streamf(j1, i1 - 1) + streamf(j1 + 1, i1 - 1) -
                         streamf(j1, i1 + 1) - streamf(j1 + 1, i1 + 1));

      advection_.set(k, globalh(j + 1, i + 1),
                     +streamf(j1 + 1, i1) - streamf(j1, i1 + 1));
    }
  advection_ *= 1. / (12. * dx * dx);
}

void StreamTransport::assembleSystem(Real dt, Real theta) {
  KmJ_ = 0.;
  KmJ_.add(1., diffusion_);
  KmJ_.add(-1., advection_);
  KmJ_.mmult(KmJmN_, bneumann_);

  mass_ = dealii::IdentityMatrix(rank_);
  mass_ /= dt;
  mass_.add(-theta, KmJmN_);
}

void StreamTransport::evolve(AthenaArray<Real>* q, Real dt, Real theta) {
  assembleSystem(dt, theta);

  for (int j = 0; j < rowsh_; ++j)
    for (int i = 0; i < colsh_; ++i) {
      int j1 = j - GhostZoneSize;
      int i1 = i - GhostZoneSize;
      int64_t k = globalh(j1, i1);
      qvec_(k) = (*q)(j, i);
    }

  KmJ_.vmult(rhs_, qvec_);
  solver_.solve(mass_, dqvec_, rhs_, dealii::PreconditionIdentity());

  for (int j = 0; j < rows_; ++j)
    for (int i = 0; i < cols_; ++i) {
      int64_t k = global(j, i);
      int j1 = j + GhostZoneSize;
      int i1 = i + GhostZoneSize;
      (*q)(j1, i1) += dqvec_(k);
    }
}
