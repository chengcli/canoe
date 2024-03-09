// athena
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/stride_iterator.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>

// exo3
#include <exo3/cubed_sphere_utility.hpp>
#include <exo3/gnomonic_equiangle.hpp>

// snap
#include "implicit_solver.hpp"

#ifdef CUBED_SPHERE
namespace cs = CubedSphereUtility;
#endif

void ImplicitSolver::SolveImplicit3D(AthenaArray<Real> &du,
                                     AthenaArray<Real> &w, Real dt) {
  auto pmb = pmy_block_;
  auto ph = pmb->phydro;
  int is, ie, js, je, ks, ke;

  // X3DIR
  if ((implicit_flag_ & (1 << 2)) && (pmb->ncells3 > 1)) {
    SetDirection(X3DIR);
    FindNeighbors();

    ks = pmb->js, js = pmb->is, is = pmb->ks;
    ke = pmb->je, je = pmb->ie, ie = pmb->ke;

    // shuffle dimension
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, i, k, j);

    // do implicit
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        if (implicit_flag_ & (1 << 3))
          FullCorrection(du_, w, dt, k, j, is, ie);
        else
          PartialCorrection(du_, w, dt, k, j, is, ie);
      }

    // shuffle back
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, i, k, j) = du_(IDN, k, j, i);
          du(IVZ, i, k, j) = du_(IVZ, k, j, i);
          du(IEN, i, k, j) = du_(IEN, k, j, i);
        }
  }

  // X2DIR
  if ((implicit_flag_ & (1 << 1)) && (pmb->ncells2 > 1)) {
    SetDirection(X2DIR);
    FindNeighbors();

    ks = pmb->is, js = pmb->ks, is = pmb->js;
    ke = pmb->ie, je = pmb->ke, ie = pmb->je;

    // shuffle dimension
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) du_(n, k, j, i) = du(n, j, i, k);

    // do implicit
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        if (implicit_flag_ & (1 << 3))
          FullCorrection(du_, w, dt, k, j, is, ie);
        else
          PartialCorrection(du_, w, dt, k, j, is, ie);
      }

    // shuffle back
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN, j, i, k) = du_(IDN, k, j, i);
          du(IVY, j, i, k) = du_(IVY, k, j, i);
          du(IEN, j, i, k) = du_(IEN, k, j, i);
        }
  }

  // X1DIR
  if (implicit_flag_ & 1) {
    SetDirection(X1DIR);
    FindNeighbors();

    ks = pmb->ks, js = pmb->js, is = pmb->is;
    ke = pmb->ke, je = pmb->je, ie = pmb->ie;

    // swap in
    du_.SwapAthenaArray(du);

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        // project velocity
#ifdef CUBED_SPHERE
        auto pco = static_cast<GnomonicEquiangle *>(pmb->pcoord);
        Real cos_theta = pco->GetCosineCell(k, j);
        Real sin_theta = pco->GetSineCell(k, j);

        for (int i = 0; i < w.GetDim1(); ++i) {
          w(IVY, k, j, i) += w(IVZ, k, j, i) * cos_theta;
          w(IVZ, k, j, i) *= sin_theta;

          // cs::CovariantToContravariant(du_.at(k,j,i), cos_theta);
          // du_(IVY, k, j, i) += du_(IVZ, k, j, i) * cos_theta;
          // du_(IVZ, k, j, i) *= sin_theta;
        }
#endif  // CUBED_SPHERE

        // do implicit
        if (implicit_flag_ & (1 << 3)) {
          FullCorrection(du_, w, dt, k, j, is, ie);
        } else {
          PartialCorrection(du_, w, dt, k, j, is, ie);
        }

        // de-project velocity
#ifdef CUBED_SPHERE
        for (int i = 0; i < w.GetDim1(); ++i) {
          w(IVZ, k, j, i) /= sin_theta;
          w(IVY, k, j, i) -= w(IVZ, k, j, i) * cos_theta;

          // du_(IVZ, k, j, i) /= sin_theta;
          // du_(IVY, k, j, i) -= du_(IVZ, k, j) * cos_theta;
          // cs::ContravariantToCovariant(du_.at(k, j, i), cos_theta);
        }
#endif  // CUBED_SPHERE
      }

    // swap out
    du_.SwapAthenaArray(du);
  }
}
