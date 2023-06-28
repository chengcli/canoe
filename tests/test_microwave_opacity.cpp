// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// snap
#include <impl.hpp>
#include <index_map.hpp>
#include <variable.hpp>

// harp
#include <harp/absorber.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// opacity
#include <opacity/Giants/microwave/mwr_absorbers.hpp>

class TestMicrowaveOpacity : public testing::Test {
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
      pmb->pindex = std::make_shared<MeshBlock::IndexMap>(pmb, pinput);
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

TEST_F(TestMicrowaveOpacity, Absorbers) {
  Radiation *prad = pmesh->my_blocks(0)->pimpl->prad;

  EXPECT_EQ(prad->GetNumBands(), 3);
  EXPECT_EQ(prad->GetBand("radio")->GetNumAbsorbers(), 3);

  Variable var;
  var.w[IDN] = 300.;
  var.w[IPR] = 1.E5;
  var.w[1] = 1.E-2;
  var.w[2] = 1.E-2;

  // NH3 absorption
  auto ab = prad->GetBand("radio")->GetAbsorber("radio-NH3");
  EXPECT_NEAR(ab->GetAttenuation(0.6, 0.6, var), 8.05853e-07, 1.E-10);

  // H2O absorption
  ab = prad->GetBand("radio")->GetAbsorber("radio-H2O");
  EXPECT_NEAR(ab->GetAttenuation(0.6, 0.6, var), 1.21548e-08, 1.E-10);

  // electron absorption
  ab = prad->GetBand("radio")->GetAbsorber("radio-Electron");
};

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  app->InstallMonitor("inversion", "main.out", "main.err");
  app->InstallMonitor("astro", "main.out", "main.err");
  app->InstallMonitor("snap", "main.out", "main.err");
  app->InstallMonitor("harp", "main.out", "main.err");

  Application::Start();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
