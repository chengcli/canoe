// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/meshblock_impl.hpp>

// harp
#include <harp/absorber.hpp>

class TestRadiation : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_radiation.inp", IOWrapper::FileMode::read);

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

TEST(TestRadiation, Construct) {
  MeshBlock *pmb = nullptr;

  Absorber ab("dummy");

  EXPECT_EQ(ab.GetName(), "dummy");
};

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  app->InstallMonitor("harp", "harp.out", "harp.err");

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
