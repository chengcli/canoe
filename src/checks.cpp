// C/C++
#include <iomanip>
#include <sstream>
#include <string>

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.h>

#include <impl.hpp>

#ifdef ENABLE_GLOG
#include <glog/logging.h>
#endif

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int k, int j, int il, int iu, int width = 1) {
  std::stringstream msg;

  msg << "rank = " << Globals::my_rank << ", " << name << " dump:" << std::endl;

  for (int i = il; i <= iu; ++i) {
    msg << "i = " << std::setw(4) << i << " |";
    msg << std::setw(4) << "(";
    for (int n = 0; n < NHYDRO; ++n) msg << var(n, k, j, i) << ", ";
    msg << ")" << std::endl;
  }

  return msg.str();
}

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int il, int iu) {
  std::stringstream msg;

  msg << "rank = " << Globals::my_rank << ", " << name << " dump:" << std::endl;

  for (int i = il; i <= iu; ++i) {
    msg << "i = " << std::setw(4) << i << " |";
    msg << std::setw(4) << "(";
    for (int n = 0; n < NHYDRO; ++n) msg << var(n, i) << ", ";
    msg << ")" << std::endl;
  }

  return msg.str();
}

void check_eos_cons2prim(AthenaArray<Real> const& prim, int k, int j, int il,
                         int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    Real const& w_d = prim(IDN, k, j, i);
    Real const& w_p = prim(IPR, k, j, i);

    int ifix = il;
    for (; ifix <= iu; ++ifix) {
      if ((prim(IDN, k, j, ifix) < 0.) || std::isnan(prim(IDN, k, j, ifix)) ||
          (prim(IPR, k, j, ifix) < 0.) || std::isnan(prim(IPR, k, j, ifix))) {
        break;
      }
    }

    LOG_IF(FATAL, std::isnan(w_d) || (w_d < 0.))
        << "ifix = " << ifix << ", "
        << print_column_table("prim-den", prim, k, j, il, iu);

    LOG_IF(FATAL, std::isnan(w_p) || (w_p < 0.))
        << "ifix = " << ifix << ", "
        << print_column_table("prim-pre", prim, k, j, il, iu);
  }
#endif  // ENABLE_GLOG
}

void check_reconstruct(AthenaArray<Real> const& wl, AthenaArray<Real> const& wr,
                       int dir, int k, int j, int il, int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    char name[80];
    snprintf(name, 80, "wl-den-%d", dir + 1);
    LOG_IF(FATAL, wl(IDN, i) < 0.) << print_column_table(name, wl, il, iu);

    snprintf(name, 80, "wr-den-%d", dir + 1);
    LOG_IF(FATAL, wr(IDN, i) < 0.) << print_column_table(name, wr, il, iu);

    snprintf(name, 80, "wl-pre-%d", dir + 1);
    LOG_IF(FATAL, wl(IPR, i) < 0.) << print_column_table(name, wl, il, iu);

    snprintf(name, 80, "wr-pre-%d", dir + 1);
    LOG_IF(FATAL, wr(IPR, i) < 0.) << print_column_table(name, wr, il, iu);
  }
#endif  // ENABLE_GLOG
}

void check_hydro_riemann_solver_flux(AthenaArray<Real> const& flux, int ivx,
                                     int k, int j, int il, int iu) {
#ifdef ENABLE_GLOG
  for (int i = il; i <= iu; ++i) {
    for (int n = 0; n < NHYDRO; ++n) {
      char name[80];
      snprintf(name, 80, "flux%d-%d", ivx, n);
      LOG_IF(FATAL, std::isnan(flux(n, k, j, i)))
          << print_column_table(name, flux, n, k, j, il, iu);
    }
  }
#endif  // ENABLE_GLOG
}

void check_decomposition(AthenaArray<Real> const& wl,
                         AthenaArray<Real> const& wr, int k, int j, int il,
                         int iu) {
  /*#ifdef ENABLE_GLOG
    for (int i = il; i <= iu; ++i) {
      LOG_IF(FATAL, wl(IDN, k, j, i) < 0.)
          << print_column_table("wl-den", wl, k, j, il, iu);

      LOG_IF(FATAL, wr(IDN, k, j, i) < 0.)
          << print_column_table("wr-den", wr, k, j, il, iu);

      LOG_IF(ERROR, wl(IPR, k, j, i) < 0.)
          << print_column_table("wl-pre", wl, k, j, il, iu);

      LOG_IF(ERROR, wr(IPR, k, j, i) < 0.)
          << print_column_table("wl-pre", wr, k, j, il, iu);
    }
  #endif  // ENABLE_GLOG*/
}

void check_implicit_cons(AthenaArray<Real> const& cons, int il, int iu, int jl,
                         int ju, int kl, int ku) {
#ifdef ENABLE_GLOG
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real const& u_d = cons(IDN, k, j, i);
        Real const& u_e = cons(IEN, k, j, i);

        LOG_IF(FATAL, std::isnan(u_d) || (u_d < 0.))
            << print_column_table("cons-den", cons, k, j, il, iu);

        LOG_IF(FATAL, std::isnan(u_e) || (u_e < 0.))
            << print_column_table("cons-eng", cons, k, j, il, iu);
      }
#endif  // ENABLE_GLOG
}

void fix_reconstruct_x1(MeshBlock* pmb, AthenaArray<Real>& wl,
                        AthenaArray<Real>& wr, AthenaArray<Real> const& w,
                        int k, int j, int il, int iu) {}

void fix_reconstruct_x2(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu) {
#ifdef ENABLE_FIX
  for (int i = il; i <= iu; ++i) {
    if (wl(IDN, i) < 0.) wl(IDN, i) = w(IDN, k, j - 1, i);
    if (wr(IDN, i) < 0.) wr(IDN, i) = w(IDN, k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = w(IPR, k, j - 1, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = w(IPR, k, j, i);
  }
#endif  // ENABLE_FIX
}

void fix_reconstruct_x3(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu) {
#ifdef ENABLE_FIX
  for (int i = il; i <= iu; ++i) {
    if (wl(IDN, i) < 0.) wl(IDN, i) = w(IDN, k - 1, j, i);
    if (wr(IDN, i) < 0.) wr(IDN, i) = w(IDN, k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = w(IPR, k - 1, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = w(IPR, k, j, i);
  }
#endif  // ENABLE_FIX
}
