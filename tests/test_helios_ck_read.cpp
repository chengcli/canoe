// external
#include <gtest/gtest.h>

#include <application/application.hpp>

// opacity
#include <opacity/absorber_ck.hpp>

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

TEST(read_helios_ck, test_case1) {}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
