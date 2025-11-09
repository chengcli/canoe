#pragma once

#include <vector>
#include <mpi.h>
#include <athena/mesh/mesh.hpp>

#include "ice_shell.hpp"

int get_mpi_rank(const MPI_Comm mpi_world = MPI_COMM_WORLD) {
  int rank;
  MPI_Comm_rank(mpi_world, &rank);
  return rank;
}


/**
 * @class SharedData
 * @brief Mimic shared memory using distributed memory
 */
template<class T>
class SharedData {
  public:
    const int n;
    const int s;
    const MPI_Comm world;
    const int root;
    const int rank;
    SharedData(const int n, const int root = 0,
        const MPI_Comm world = MPI_COMM_WORLD);

    inline T get(int i) const {
      return value[i];
    }

    inline void set(int i, T x) {
      value[i] = x;
      modified[i] = true;
    }

    void share(void);
  private:
    std::vector<T> value;
    std::vector<int> modified;
    std::vector<T> global_value;
    std::vector<int> global_modified;

    inline static int mpi_world_size(const MPI_Comm mpi_world) {
      int world_size;
      MPI_Comm_size(mpi_world, &world_size);
      return world_size;
    }

    inline static int mpi_rank(const MPI_Comm mpi_world) {
      int rank;
      MPI_Comm_rank(mpi_world, &rank);
      return rank;
    }

    inline static int global_allocation_size(
        const int n, const int root, const MPI_Comm mpi_world) {
      return (mpi_rank(mpi_world) == root) ? n * mpi_world_size(mpi_world) : 0;
    }
};

template<class T>
SharedData<T>::SharedData(const int n, const int root, const MPI_Comm world):
  n(n), world(world), s(mpi_world_size(world)),
  rank(mpi_rank(world)), root(root),
  value(n), modified(n, false),
  global_value(global_allocation_size(n, root, world)),
  global_modified(global_allocation_size(n, root, world)) {
}

inline auto mpi_element_type(std::vector<double> x) {
  return MPI_DOUBLE;
}

inline auto mpi_element_type(std::vector<int> x) {
  return MPI_INT;
}

template<class T>
void SharedData<T>::share(void) {

  MPI_Gather(value.data(), n, mpi_element_type(value),
      global_value.data(), n, mpi_element_type(value),
      root, MPI_COMM_WORLD);

  MPI_Gather(modified.data(), n, mpi_element_type(modified),
      global_modified.data(), n, mpi_element_type(modified),
      root, MPI_COMM_WORLD);

  if (rank == root) {
    for (int j = 0, k = 0; j < s; ++j) {
      for (int i = 0; i < n; ++i, ++k) {
        if (global_modified[k]) {
            value[i] = global_value[k];
        }
      }
    }
  }

  MPI_Bcast(value.data(), n, mpi_element_type(value), root, MPI_COMM_WORLD);
  std::fill(modified.begin(), modified.end(), false);
}

template<typename Real>
auto init_ice_boundary_model (
    MeshBlock * pmb, const Real ice_max_x1, const Real ice_min_x2) {
  const Real vapor_p3 = 611.7; // Pascal
  const Real vapor_t3 = 273.16; // Kelvin
  const Real vapor_gas_const = 461.5; // J kg-1 K-1
  const Real vapor_adiabatic_index = 1.4;
  const Real vapor_beta = 24.845;
  const Real vapor_delta = 4.986;

  const Real outer_surface_temp_eff = 67.; // Kelvin
  const Real stefan_boltzmann_const = 5.67e-8;

  const Real ocean_temp = 273.16; // Kelvin
  const Real ice_k = 651.; // Watts per meter

  IceShell::VaporCondensation<Real> cond
    (vapor_p3, vapor_t3, vapor_gas_const, vapor_adiabatic_index,
       vapor_beta, vapor_delta);
  IceShell::OuterSurfaceRadiation<Real> rad
    (outer_surface_temp_eff, stefan_boltzmann_const);
  IceShell::TemperatureDependentDiffusivity<Real> tdd
    (ocean_temp, ice_k);

  int lx1 = pmb->loc.lx1;
  int lx2 = pmb->loc.lx2;
  int nx1 = pmb->block_size.nx1;
  int nx2 = pmb->block_size.nx2;
  int g_nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int g_nx2 = pmb->pmy_mesh->mesh_size.nx2;

  SharedData<Real> mesh_dx1f(g_nx1);
  SharedData<Real> mesh_dx2f(g_nx2);
  SharedData<Real> mesh_x1v(g_nx1);
  SharedData<Real> mesh_x2v(g_nx2);

  for (int i = pmb->is, ig=lx1*nx1; i <= pmb->ie; ++ig, ++i) {
    mesh_dx1f.set(ig, pmb->pcoord->dx1f(i));
    mesh_x1v.set(ig, pmb->pcoord->x1v(i));
  }

  for (int j = pmb->js, jg=lx2*nx2; j <= pmb->je; ++jg, ++j) {
    mesh_dx2f.set(jg, pmb->pcoord->dx2f(j));
    mesh_x2v.set(jg, pmb->pcoord->x2v(j));
  }

  mesh_dx1f.share();
  mesh_dx2f.share();
  mesh_x1v.share();
  mesh_x2v.share();

  std::vector<Real> dx;
  std::vector<Real> dz;

  for (int ig = 0; ig < mesh_dx1f.n; ++ig) {
    if (mesh_x1v.get(ig) < ice_max_x1) {
      dz.push_back(mesh_dx1f.get(ig));
    }
  }

  for (int jg = 0; jg < mesh_dx2f.n; ++jg) {
    if (mesh_x2v.get(jg) > ice_min_x2) {
      dx.push_back(mesh_dx2f.get(jg));
    }
  }

  IceShell::IceBoundaryModel<Real> ice_boundary_model
    (dx, dz, cond, rad, tdd);

  if (get_mpi_rank() == 0) {
    std::cout << "Ice Model: nx = " << dx.size()
      << "; nz = " << dz.size() << std::endl;
  }

  return ice_boundary_model;
}

template<typename R1, typename R2, typename R3>
inline bool fclose(R1 x, R2 x0, R3 abs_tol) {
  return std::abs(x - x0) < abs_tol;
}

template<typename Real>
inline bool meshblock_is_right_ice(MeshBlock *pmb,
    Real ice_max_x1, Real ice_min_x2) {
  return (
    fclose(pmb->block_size.x2max, ice_min_x2, 1e-6)
    && pmb->block_size.x1max < (ice_max_x1 + 1e-6)
  );
}

template<typename Real>
inline bool meshblock_is_bottom_ice(MeshBlock *pmb,
    Real ice_max_x1, Real ice_min_x2) {
  return (
    fclose(pmb->block_size.x1min, ice_max_x1, 1e-6)
    && pmb->block_size.x2min > (ice_min_x2 - 1e-6)
  );
}

template<typename Real>
class AirIceCoupler {
  public:
    AirIceCoupler(MeshBlock *pmb,
      const Real ice_max_x1, const Real ice_min_x2, const Real drag_coef,
      const int i_vapor,
      const Real abs_tol, const Real max_iter):
        i_vapor(i_vapor),
        drag_coef(drag_coef),
        abs_tol(abs_tol),
        max_iter(max_iter),
        is_root(get_mpi_rank() == 0),
        is_right_ice(meshblock_is_right_ice(pmb, ice_max_x1, ice_min_x2)),
        is_bottom_ice(meshblock_is_bottom_ice(pmb, ice_max_x1, ice_min_x2)),
        ibm(init_ice_boundary_model(pmb, ice_max_x1, ice_min_x2)),
        i_offset(pmb->loc.lx1 * pmb->block_size.nx1 - pmb->is),
        j_offset(pmb->loc.lx2 * pmb->block_size.nx2 - pmb->js
              - pmb->pmy_mesh->mesh_size.nx2 + ibm.nx),
        air_t(ibm.nx, ibm.nz),
        vapor_p(ibm.nx, ibm.nz),
        ice_t(ibm.nx, ibm.nz),
        air_t_side(ibm.nz),
        air_t_top(ibm.nx),
        vapor_p_side(ibm.nz),
        vapor_p_top(ibm.nx),
        ice_t_side(ibm.nz),
        ice_t_top(ibm.nx) {
    }

    void solve(MeshBlock *pmb, AthenaArray<Real> const &w);

    inline int ice_i(int i) const {
      return i + i_offset;
    }

    inline int ice_j(int j) const {
      return j + j_offset;
    }

    inline auto ice_forcing_on_air(MeshBlock *pmb,
          int k, int j, int i) const {

      struct AirTendency {
        Real H2O;
        Real rho_u1;
        Real rho_u2;
        Real en;
      } g = {0, 0, 0, 0};

      const Real rho = pmb->phydro->w(IDN, k, j, i);
      const Real u1 = pmb->phydro->w(IVX, k, j, i);
      const Real u2 = pmb->phydro->w(IVY, k, j, i);

      // ice to the right
      if (is_right_ice && j == pmb->je) {

        // condensation and evaporation
        const int l = ice_i(i);
        g.H2O -= (
          condensation_rate(
            air_t_side.get(l), vapor_p_side.get(l), ice_t_side.get(l)
          ) / pmb->pcoord->dx2f(j)
        );

        // friction drag
        g.rho_u1 -= (
          2 * drag_coef * rho * std::abs(u1) * u1
          / pmb->pcoord->dx2f(j)
        );
      }

      // ice at the bottom
      if (is_bottom_ice && i == pmb->is) {

        // condensation and evaporation
        const int l = ice_j(j);
        g.H2O -= (
          condensation_rate(
            air_t_top.get(l), vapor_p_top.get(l), ice_t_top.get(l)
          ) / pmb->pcoord->dx1f(i)
        );

        // friction drag
        g.rho_u2 -= (
          2 * drag_coef * rho * std::abs(u2) * u2
          / pmb->pcoord->dx1f(i)
        );
      }

      // momentom and energy exchange associated with mass exchange
      const auto pthermo = Thermodynamics::GetInstance();
      const Real ie = (
        (pthermo->GetRd() / (pthermo->GetGammad() - 1.))
        * pthermo->GetCvRatio(i_vapor)
        * pthermo->GetTemp(pmb->phydro->w.at(k, j, i))
      );

      if (g.H2O < 0) {
        g.rho_u1 += g.H2O * u1;
        g.rho_u2 += g.H2O * u2;
        g.en += g.H2O * (0.5 * (u1 * u1 + u2 * u2) + ie);
      } else {
        g.en += g.H2O * ie;
      }

      return g;
    }

    inline Real adjacent_ice_t(MeshBlock *pmb, int i, int j) {
      Real g = -1.;
      if (is_right_ice && j == pmb->je) {
        int l = ice_i(i);
        g = ice_t_side.get(l);
      }
      if (is_bottom_ice && i == pmb->is) {
        int l = ice_j(j);
        g = ice_t_top.get(l);
      }
      return g;
    }

  private:
    const int i_vapor;
    const Real drag_coef;
    const Real abs_tol;
    const int max_iter;
    const bool is_root;
    const bool is_right_ice;
    const bool is_bottom_ice;
    IceShell::IceBoundaryModel<Real> ibm;
    const int i_offset;
    const int j_offset;
    IceShell::BoundaryValue<Real> air_t;
    IceShell::BoundaryValue<Real> vapor_p;
    IceShell::BoundaryValue<Real> ice_t;
    SharedData<Real> air_t_side;
    SharedData<Real> air_t_top;
    SharedData<Real> vapor_p_side;
    SharedData<Real> vapor_p_top;
    SharedData<Real> ice_t_side;
    SharedData<Real> ice_t_top;

    inline Real condensation_rate(
        Real air_t, Real vapor_p, Real ice_t) const {
      return ibm.ice_air_boundary.cond.net_vapor_flux(ice_t, air_t, vapor_p);
    }
};

template<typename Real>
Real get_air_t(AthenaArray<Real> const &w, int k, int j, int i) {
  auto pthermo = Thermodynamics::GetInstance();
  return pthermo->GetTemp(w.at(k, j, i));
}

template<typename Real>
Real get_vapor_p(AthenaArray<Real> const &w, int k, int j, int i, int iv) {
  auto pthermo = Thermodynamics::GetInstance();
  return (
    w(IDN, k, j, i) * w(iv, k, j, i)
    * pthermo->GetTemp(w.at(k, j, i))
    * pthermo->GetRd() * pthermo->GetInvMuRatio(iv)
  );
}

template<typename Real>
void AirIceCoupler<Real>::solve(MeshBlock *pmb, AthenaArray<Real> const &w) {
  int k = pmb->ks; // two-dimensional flow
  if (is_right_ice) {
    int j = pmb->je;
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      int l = ice_i(i);
      air_t_side.set(l, get_air_t(w, k, j, i));
      vapor_p_side.set(l, get_vapor_p(w, k, j, i, i_vapor));
    }
  }
  if (is_bottom_ice) {
    int i = pmb->is;
    for (int j = pmb->js; j <= pmb->je; ++j) {
      int l = ice_j(j);
      air_t_top.set(l, get_air_t(w, k, j, i));
      vapor_p_top.set(l, get_vapor_p(w, k, j, i, i_vapor));
    }
  }
  air_t_side.share();
  air_t_top.share();
  vapor_p_side.share();
  vapor_p_top.share();
  if (is_root) {
    for (int i = 0; i < ibm.nz; ++i) {
      air_t.set_side(i, air_t_side.get(i));
      vapor_p.set_side(i, vapor_p_side.get(i));
    }
    for (int i = 0; i < ibm.nx; ++i) {
      air_t.set_top(i, air_t_top.get(i));
      vapor_p.set_top(i, vapor_p_top.get(i));
    }
    ibm.solve(ice_t, air_t, vapor_p, abs_tol, max_iter);
    for (int i = 0; i < ibm.nz; ++i) {
      ice_t_side.set(i, ice_t.get_side(i));
    }
    for (int i = 0; i < ibm.nx; ++i) {
      ice_t_top.set(i, ice_t.get_top(i));
    }
  }
  ice_t_side.share();
  ice_t_top.share();
}
