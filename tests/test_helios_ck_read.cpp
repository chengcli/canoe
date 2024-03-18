// C/C++
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// external
#include <gtest/gtest.h>

#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// opacity
#include <opacity/absorber_ck.hpp>

// canoe
#include <impl.hpp>

std::string data_folder = "ck_data_01242024/ck/";

class TestHelioCK : public testing::Test {
 protected:
  Mesh* pmesh;
  ParameterInput* pinput;

  virtual void SetUp() {
    IOWrapper infile;
    infile.Open("test_helios_ck.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    Thermodynamics::InitFromAthenaInput(pinput);

    // set up mesh
    int restart = false;
    int mesh_only = false;
    pmesh = new Mesh(pinput, mesh_only);

    // set up components
    for (int b = 0; b < pmesh->nblocal; ++b) {
      MeshBlock* pmb = pmesh->my_blocks(b);
      pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
    }
  }

  virtual void TearDown() {
    Thermodynamics::Destroy();
    IndexMap::Destroy();
    delete pinput;
    delete pmesh;
  }
};

TEST(LoadCoefficient, bid_is_0) {
  auto app = Application::GetInstance();
  std::string fname = "PM_ck_HELIOSK_cond_11_nOPT_wcia.txt";
  auto file = app->FindResource(data_folder + fname);
  HeliosCK PM_ck_HELIOSK_cond_11_nOPT_wcia("PM_ck_HELIOSK_cond_11_nOPT_wcia");
  PM_ck_HELIOSK_cond_11_nOPT_wcia.LoadCoefficient(file, 0);
  EXPECT_EQ(1, 1);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
