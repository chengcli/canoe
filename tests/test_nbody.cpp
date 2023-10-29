// C/C++
#include <iostream>

// external
#include <gtest/gtest.h>
#include <mpi.h>
#include <omp.h>

// application
#include <application/application.hpp>

// athena
#include <athena/globals.hpp>

// pvfmm
#include <pvfmm/pvfmm.hpp>

#include "pvfmm_utils.hpp"

using vecd = std::vector<double>;
using veci = std::vector<int>;

void nbody_direct(vecd& sl_coord, vecd& sl_den, vecd& dl_coord, vecd& dl_den,
                  vecd& trg_coord, vecd& trg_value,
                  const pvfmm::Kernel<double>& kernel_fn) {
  long long n_sl = sl_coord.size() / PVFMM_COORD_DIM;
  long long n_dl = dl_coord.size() / PVFMM_COORD_DIM;
  long long n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  long long n_trg_glb = 0;

  MPI_Allreduce(&n_trg, &n_trg_glb, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  // ker_dim is a MxN matrix
  // For Laplace kernel, M = N = 1
  vecd glb_trg_coord(n_trg_glb * PVFMM_COORD_DIM);
  vecd glb_trg_value(n_trg_glb * kernel_fn.ker_dim[1], 0);

  veci recv_disp(Globals::nranks, 0);

  // gather trg_coord and trg_value to all ranks
  int send_cnt = n_trg * PVFMM_COORD_DIM;
  veci recv_cnts(Globals::nranks);
  MPI_Allgather(&send_cnt, 1, MPI_INT, recv_cnts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);
  MPI_Allgatherv(trg_coord.data(), send_cnt, MPI_DOUBLE, glb_trg_coord.data(),
                 recv_cnts.data(), recv_disp.data(), MPI_DOUBLE,
                 MPI_COMM_WORLD);

  // evaluate target potential
  vecd glb_trg_value_(n_trg_glb * kernel_fn.ker_dim[1], 0);
  int nthreads = omp_get_max_threads();

#pragma omp parallel
  for (int i = 0; i < nthreads; ++i) {
    size_t a = (i * n_trg_glb) / nthreads;
    size_t b = ((i + 1) * n_trg_glb) / nthreads;

    if (kernel_fn.ker_poten != NULL) {
      kernel_fn.ker_poten(sl_coord.data(), n_sl, sl_den.data(), 1,
                          glb_trg_coord.data() + a * PVFMM_COORD_DIM, b - a,
                          glb_trg_value_.data() + a * kernel_fn.ker_dim[1],
                          NULL);
    }

    if (kernel_fn.dbl_layer_poten != NULL) {
      kernel_fn.dbl_layer_poten(
          dl_coord.data(), n_dl, dl_den.data(), 1,
          glb_trg_coord.data() + a * PVFMM_COORD_DIM, b - a,
          glb_trg_value_.data() + a * kernel_fn.ker_dim[1], NULL);
    }

    // gather contributions from all ranks
    MPI_Allreduce(glb_trg_value_.data(), glb_trg_value.data(),
                  glb_trg_value.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // set local target values
    int disp = recv_disp[Globals::my_rank] / PVFMM_COORD_DIM;
    trg_value.assign(
        glb_trg_value.begin() + disp * kernel_fn.ker_dim[1],
        glb_trg_value.begin() + (disp + n_trg) * kernel_fn.ker_dim[1]);
  }
}

TEST(fmm_laplace, potential) {
  // set kernel
  const pvfmm::Kernel<double>& kernel_fn =
      pvfmm::LaplaceKernel<double>::potential();
  auto comm = MPI_COMM_WORLD;
  int multipole_order = 10;

  // Create target and source vectors
  vecd trg_coord = {0.1, 0, 0, 0.2, 0, 0, 0.3, 0, 0};
  vecd sl_coord = {0, 0, 0};
  vecd dl_coord;

  // Dimension of the problem
  size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  size_t n_sl = sl_coord.size() / PVFMM_COORD_DIM;
  size_t n_dl = dl_coord.size() / PVFMM_COORD_DIM;

  // Set source charges
  vecd sl_den = {1.};
  vecd dl_den;

  // Create memory-manager (optional)
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  {  // FMM N-Body
    // Construct tree
    size_t max_pts = 600;
    auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den, trg_coord,
                                  MPI_COMM_WORLD, max_pts, pvfmm::FreeSpace);

    // Load matrices
    pvfmm::PtFMM<double> matrices(&mem_mgr);
    matrices.Initialize(multipole_order, MPI_COMM_WORLD, &kernel_fn);

    // FMM Setup
    tree->SetupFMM(&matrices);

    // Run FMM
    vecd trg_value;
    PtFMM_Evaluate(tree, trg_value, n_trg);

    EXPECT_NEAR(trg_value[0], 1. / (4. * M_PI * 0.1), 1e-10);
    EXPECT_NEAR(trg_value[1], 1. / (4. * M_PI * 0.2), 1e-10);
    EXPECT_NEAR(trg_value[2], 1. / (4. * M_PI * 0.3), 1e-10);

    // Free memory
    delete tree;
  }

  {  // Direct N-Body
    vecd trg_value_direct(n_trg * kernel_fn.ker_dim[1]);
    nbody_direct(sl_coord, sl_den, dl_coord, dl_den, trg_coord,
                 trg_value_direct, kernel_fn);
    EXPECT_NEAR(trg_value_direct[0], 1. / (4. * M_PI * 0.1), 1e-10);
  }
}

TEST(fmm_laplace, example1) {
  // Set kernel.
  const pvfmm::Kernel<double>& kernel_fn =
      pvfmm::LaplaceKernel<double>::gradient();
  auto comm = MPI_COMM_WORLD;
  int N = 10;
  int mult_order = 10;

  // Create target and source vectors.
  vecd trg_coord = point_distrib<double>(RandUnif, N, comm);
  vecd sl_coord = point_distrib<double>(RandUnif, N, comm);
  vecd dl_coord = point_distrib<double>(RandUnif, 0, comm);

  size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  size_t n_sl = sl_coord.size() / PVFMM_COORD_DIM;
  size_t n_dl = dl_coord.size() / PVFMM_COORD_DIM;

  // Set source charges.
  vecd sl_den(n_sl * kernel_fn.ker_dim[0]);
  vecd dl_den(n_dl * (kernel_fn.ker_dim[0] + PVFMM_COORD_DIM));
  for (size_t i = 0; i < sl_den.size(); i++) sl_den[i] = drand48();
  for (size_t i = 0; i < dl_den.size(); i++) dl_den[i] = drand48();

  // Create memory-manager (optional)
  pvfmm::mem::MemoryManager mem_mgr(10000000);

  // Construct tree.
  size_t max_pts = 600;
  auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den, trg_coord,
                                comm, max_pts, pvfmm::FreeSpace);

  // Load matrices.
  pvfmm::PtFMM<double> matrices(&mem_mgr);
  matrices.Initialize(mult_order, comm, &kernel_fn);

  // FMM Setup
  tree->SetupFMM(&matrices);

  // Run FMM
  vecd trg_value;
  PtFMM_Evaluate(tree, trg_value, n_trg);

  // Re-run FMM
  tree->ClearFMMData();
  for (size_t i = 0; i < sl_den.size(); i++) sl_den[i] = drand48();
  for (size_t i = 0; i < dl_den.size(); i++) dl_den[i] = drand48();
  PtFMM_Evaluate(tree, trg_value, n_trg, &sl_den, &dl_den);
}

int main(int argc, char** argv) {
  Application::Start(argc, argv);
  omp_set_num_threads(1);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();
  return result;
}
