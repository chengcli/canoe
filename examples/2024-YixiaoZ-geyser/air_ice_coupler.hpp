#pragma once

#include <vector>
#include <mpi.h>
#include <athena/mesh/mesh.hpp>

#include "ice_shell.hpp"

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

template<class T>
void SharedData<T>::share(void) {

  MPI_Gather(value.data(), n, MPI_DOUBLE,
      global_value.data(), n, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Gather(modified.data(), n, MPI_INT,
      global_modified.data(), n, MPI_INT, root, MPI_COMM_WORLD);

  if (rank == root) {
    for (int j = 0, k = 0; j < s; ++j) {
      for (int i = 0; i < n; ++i, ++k) {
        if (global_modified[k]) {
            value[i] = global_value[k];
        }
      }
    }
  }

  MPI_Bcast(value.data(), n, MPI_DOUBLE, root, MPI_COMM_WORLD);
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

  std::cout << "Ice Model: nx = " << dx.size()
    << "; nz = " << dz.size() << std::endl;

  return ice_boundary_model;
}

template<typename R1, typename R2, typename R3>
inline bool fclose(R1 x, R2 x0, R3 abs_tol) {
  return std::abs(x - x0) < abs_tol;
}

template<typename Real>
class AirIceCoupler {
  public:
    AirIceCoupler(MeshBlock *pmb,
      const Real ice_max_x1, const Real ice_min_x2):
        is_right_ice(fclose(pmb->block_size.x2max, ice_min_x2, 1e-6)),
        is_bottom_ice(fclose(pmb->block_size.x1min, ice_max_x1, 1e-6)),
        g_i(pmb->loc.lx1 * pmb->block_size.nx1),
        g_j(pmb->loc.lx2 * pmb->block_size.nx2),
        ibm(init_ice_boundary_model(pmb, ice_max_x1, ice_min_x2)),
        air_t(ibm.nx, ibm.nz),
        ice_t(ibm.nx, ibm.nz),
        air_t_side(ibm.nz),
        air_t_top(ibm.nx),
        ice_t_side(ibm.nz),
        ice_t_top(ibm.nx) {
    }
  private:
    const bool is_right_ice;
    const bool is_bottom_ice;
    const int g_i;
    const int g_j;
    IceShell::IceBoundaryModel<Real> ibm;
    IceShell::BoundaryValue<Real> air_t;
    IceShell::BoundaryValue<Real> ice_t;
    SharedData<Real> air_t_side;
    SharedData<Real> air_t_top;
    SharedData<Real> ice_t_side;
    SharedData<Real> ice_t_top;
};
