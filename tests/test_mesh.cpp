// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// canoe
#include <impl.hpp>

class TestMesh : public testing::Test {
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
nx1         = 4         # Number of zones in X1-direction
x1min       = -0.5      # minimum value of X1
x1max       = 0.5       # maximum value of X1
ix1_bc      = outflow   # Inner-X1 boundary condition flag
ox1_bc      = outflow   # Outer-X1 boundary condition flag

nx2         = 4         # Number of zones in X2-direction
x2min       = -0.5      # minimum value of X2
x2max       = 0.5       # maximum value of X2
ix2_bc      = periodic  # Inner-X2 boundary condition flag
ox2_bc      = periodic  # Outer-X2 boundary condition flag

nx3         = 1         # Number of zones in X3-direction
x3min       = -0.5      # minimum value of X3
x3max       = 0.5       # maximum value of X3
ix3_bc      = periodic  # Inner-X3 boundary condition flag
ox3_bc      = periodic  # Outer-X3 boundary condition flag

<hydro>
gamma       = 1.4

<problem>
dens  = 1.0
pres  = 1.E5
qH2O  = 0.1
qNH3  = 0.01
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

TEST_F(TestMesh, Problem) {
  // read in parameters
  Real dens = pinput->GetReal("problem", "dens");
  Real pres = pinput->GetReal("problem", "pres");
  Real qH2O = pinput->GetReal("problem", "qH2O");
  Real qNH3 = pinput->GetReal("problem", "qNH3");

  EXPECT_EQ(dens, 1.0);
  EXPECT_EQ(pres, 1.0e5);
  EXPECT_EQ(qH2O, 0.1);
  EXPECT_EQ(qNH3, 0.01);
};

TEST_F(TestMesh, Mesh) {
  EXPECT_EQ(pmesh->mesh_size.nx1, 4);
  EXPECT_EQ(pmesh->mesh_size.x1min, -0.5);
  EXPECT_EQ(pmesh->mesh_size.x1max, 0.5);

  EXPECT_EQ(pmesh->mesh_size.nx2, 4);
  EXPECT_EQ(pmesh->mesh_size.x2min, -0.5);
  EXPECT_EQ(pmesh->mesh_size.x2max, 0.5);

  EXPECT_EQ(pmesh->mesh_size.nx3, 1);
  EXPECT_EQ(pmesh->mesh_size.x3min, -0.5);
  EXPECT_EQ(pmesh->mesh_size.x3max, 0.5);
};

int main(int argc, char **argv) {
  Application::Start(argc, argv);
  testing::InitGoogleTest(&argc, argv);

  auto app = Application::GetInstance();

  app->InstallMonitor("inversion", "main.out", "main.err");
  app->InstallMonitor("astro", "main.out", "main.err");
  app->InstallMonitor("snap", "main.out", "main.err");
  app->InstallMonitor("harp", "main.out", "main.err");

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
