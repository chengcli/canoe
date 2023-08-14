#ifndef CHEMISTRY_SOLVER_HPP
#define CHEMISTRY_SOLVER_HPP

// #include "../math/linalg.h"
//  Athena++ headers
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"
#include "../math/eigen335/Eigen/LU"

template <int N>
class ChemistrySolver {
 public:
  // needed for Eigen small matrix
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // data
  enum { Size = N };

  // functions
  ChemistrySolver() { I_.setIdentity(); }

  template <typename T1, typename T2>
  T1 solveBDF1(T1 const& Rate, T2 const& Jac, Real dt) {
    A_ = I_ / dt - Jac;
    if (N <= 4)
      return A_.inverse() * Rate;
    else
      return A_.partialPivLu().solve(Rate);
  }

  template <typename T1, typename T2>
  T1 solveTRBDF2(T1 const& Rate, T2 const& Jac, Real dt) {
    int gamma = 2. - sqrt(2.);
    A_ = I_ - gamma / 2. * dt * Jac;
    B_ = dt * (1. - gamma / 2.) * Rate + dt * gamma / 2. * A_ * Rate;
    A_ = A_ * A_;
    if (N <= 4)
      return A_.inverse() * B_;
    else
      return A_.partialPivLu().solve(B_);
  }

  template <typename T1, typename T2>
  T1 solveTRBDF2Blend(Real c[], T1 const& Rate, T2 const& Jac, int const indx[],
                      Real dt) {
    // BDF1 solver
    Eigen::Matrix<Real, N, 1> sol = BDF1(Rate, Jac, dt);
    for (int n = 0; n < Size; ++n) S1_(n) = c[indx[n]] + sol(n);

    // TR-BDF2 solver
    sol = TRBDF2(Rate, Jac, dt);
    for (int n = 0; n < Size; ++n) S2_(n) = c[indx[n]] + sol(n);

    // Blend solutions
    Real alpha = 1.;
    for (int n = 0; n < Size; ++n)
      if (S2_(n) < 0.) alpha = std::min(alpha, S1_(n) / (S1_(n) - S2_(n)));
    for (int n = 0; n < Size; ++n) {
      sol(n) = (1. - alpha) * S1_(n) + alpha * S2_(n) - c[indx[n]];
      c[indx[n]] += sol(n);
    }
    return sol;
  }

 private:
  // scratch array
  Eigen::Matrix<Real, N, N> A_, I_;
  Eigen::Matrix<Real, N, 1> B_, S1_, S2_;
};

#endif
