// C/C++
#include <string>

// athena
#include <athena/athena.hpp>
#include <athena/globals.hpp>

class MeshBlock;

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int k, int j, int il, int iu, int width = 1);

std::string print_column_table(std::string name, AthenaArray<Real> const& var,
                               int il, int iu);

void check_eos_cons2prim(AthenaArray<Real> const& prim, int k, int j, int il,
                         int iu);

void fix_eos_cons2prim(MeshBlock* pmb, AthenaArray<Real>& prim,
                       AthenaArray<Real>& cons, int k, int j, int il, int iu);

void check_reconstruct(AthenaArray<Real> const& wl, AthenaArray<Real> const& wr,
                       int ivx, int k, int j, int il, int iu);

void fix_reconstruct_x1(MeshBlock* pmb, AthenaArray<Real>& wl,
                        AthenaArray<Real>& wr, AthenaArray<Real> const& w,
                        int k, int j, int il, int iu);

void fix_reconstruct_x2(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu);

void fix_reconstruct_x3(AthenaArray<Real>& wl, AthenaArray<Real>& wr,
                        AthenaArray<Real> const& w, int k, int j, int il,
                        int iu);

void check_hydro_riemann_solver_flux(AthenaArray<Real> const& flux, int ivx,
                                     int k, int j, int il, int iu);

void check_decomposition(AthenaArray<Real> const& wl,
                         AthenaArray<Real> const& wr, int k, int j, int il,
                         int iu);

void check_implicit_cons(AthenaArray<Real> const& cons, int il, int iu, int jl,
                         int ju, int kl, int ku);

void fix_implicit_cons(MeshBlock* pmb, AthenaArray<Real>& cons, int il, int iu,
                       int jl, int ju, int kl, int ku);
