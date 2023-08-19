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

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

class TestMicrophysics : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;

  virtual void SetUp() {
    IOWrapper infile;
    infile.Open("test_microphysics.inp", IOWrapper::FileMode::read);

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
    delete pinput;
    delete pmesh;
    Thermodynamics::Destroy();
    IndexMap::Destroy();
  }
};

TEST_F(TestMicrophysics, initialize) {}

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
