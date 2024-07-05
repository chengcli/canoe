// external
#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/reconstruct/reconstruction.hpp>

// torch
#include <torch/torch.h>

// canoe
#include <impl.hpp>
#include <snap/athena_arrays.hpp>

class TestReconstruct : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;
  char fname[80] = "/tmp/tempfile.XXXXXX";

  void CreateInputFile() {
    const char *mesh_config = R"(
<time>
cfl_number  = 1
tlim        = 0.1

<mesh>
nx1         = 100         # Number of zones in X1-direction
x1min       = -0.5      # minimum value of X1
x1max       = 0.5       # maximum value of X1
ix1_bc      = outflow   # Inner-X1 boundary condition flag
ox1_bc      = outflow   # Outer-X1 boundary condition flag

nx2         = 1000         # Number of zones in X2-direction
x2min       = -0.5      # minimum value of X2
x2max       = 0.5       # maximum value of X2
ix2_bc      = periodic  # Inner-X2 boundary condition flag
ox2_bc      = periodic  # Outer-X2 boundary condition flag

nx3         = 1000         # Number of zones in X3-direction
x3min       = -0.5      # minimum value of X3
x3max       = 0.5       # maximum value of X3
ix3_bc      = periodic  # Inner-X3 boundary condition flag
ox3_bc      = periodic  # Outer-X3 boundary condition flag

<hydro>
gamma       = 1.4
)";
    // write to file
    mkstemp(fname);
    std::ofstream outfile(fname);
    outfile << mesh_config;
  }

  virtual void SetUp() {
    CreateInputFile();

    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open(fname, IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    // set up mesh
    int restart = false;
    int mesh_only = false;
    pmesh = new Mesh(pinput, mesh_only);

    // set up components
    for (int b = 0; b < pmesh->nblocal; ++b) {
      MeshBlock *pmb = pmesh->my_blocks(b);
      pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
    }

    pmesh->Initialize(restart, pinput);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    delete pinput;
    delete pmesh;
    std::remove(fname);
  }
};

/*TEST_F(TestReconstruct, test1) {
  auto pmb = pmesh->my_blocks(0);
  int nc1 = pmb->ncells1;
  int nc2 = pmb->ncells2;
  int nc3 = pmb->ncells3;

  AthenaArray<float> w(NHYDRO, nc3, nc2, nc1);
  w.toDevice(torch::kMPS);
  w.tensor().normal_(0, 1);
  w.fromDevice();

  for (int n = 0; n < NHYDRO; ++n) {
    for (int k = 0; k < nc3; ++k) {
      std::cout << "n = " << n << ", k = " << k << std::endl;
      std::cout << "---------" << std::endl;
      for (int j = 0; j < nc2; ++j) {
        for (int i = 0; i < nc1; ++i) {
          std::cout << w(n, k, j, i) << ", ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl << "---------" << std::endl;
    }
    std::cout << std::endl;
  }
}*/

TEST_F(TestReconstruct, test_x1) {
  auto pmb = pmesh->my_blocks(0);

  AthenaArray<float> wl3d, wr3d, w;

  int nc1 = pmb->ncells1;
  int nc2 = pmb->ncells2;
  int nc3 = pmb->ncells3;

  w.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
  wl3d.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
  wr3d.NewAthenaArray(NHYDRO, nc3, nc2, nc1);

  w.toDevice(torch::kMPS);
  wl3d.toDevice(torch::kMPS);
  wr3d.toDevice(torch::kMPS);

  w.tensor().normal_(0, 1);

  auto start = std::chrono::high_resolution_clock::now();
  if (NGHOST > 2) {
    pmb->precon->Weno5X1(pmb->is - 1, pmb->ie + 1, w.tensor(), wl3d.tensor(),
                         wr3d.tensor());
    pmb->precon->Weno5X2(pmb->js - 1, pmb->je + 1, w.tensor(), wl3d.tensor(),
                         wr3d.tensor());
    pmb->precon->Weno5X3(pmb->ks - 1, pmb->ke + 1, w.tensor(), wl3d.tensor(),
                         wr3d.tensor());
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Time taken by test body: " << elapsed.count() << " seconds"
            << std::endl;

  // wl3d.fromDevice();
  // wr3d.fromDevice();
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
