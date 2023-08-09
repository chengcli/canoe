// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// harp
#include <harp/absorber.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

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

    IndexMap::InitFromAthenaInput(pinput);
    Thermodynamics::InitFromAthenaInput(pinput);

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
    Thermodynamics::Destroy();
    IndexMap::Destroy();
  }
};

TEST_F(TestRadiation, Species) {
  auto pindex = IndexMap::GetInstance();
  EXPECT_EQ(pindex->GetSpeciesId("vapor.H2O"), 1);
  EXPECT_EQ(pindex->GetSpeciesId("vapor.NH3"), 2);
}

TEST_F(TestRadiation, Radiation) {
  auto prad = pmesh->my_blocks(0)->pimpl->prad;

  EXPECT_EQ(prad->GetNumBands(), 3);
  EXPECT_EQ(prad->GetBand(0)->GetNumAbsorbers(), 6);
  EXPECT_EQ(prad->GetBand(0)->GetAbsorber(0)->GetName(), "H2-H2-CIA");
  EXPECT_EQ(prad->GetBand(0)->GetAbsorber(1)->GetName(), "H2-He-CIA");
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
