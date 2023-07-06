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

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_mesh.inp", IOWrapper::FileMode::read);

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

  EXPECT_EQ(pmesh->mesh_size.nx2, 2);
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
