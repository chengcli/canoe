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

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

class TestThermodynamics : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_thermodynamics.inp", IOWrapper::FileMode::read);

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
    Thermodynamics::InitFromAthenaInput(pinput);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    delete pinput;
    delete pmesh;
  }
};

TEST_F(TestThermodynamics, dry) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetRd(), 3777., 1e-8);
  EXPECT_NEAR(pthermo->GetMu(0), 0.00220134, 1e-8);
}

TEST_F(TestThermodynamics, water_vapor) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetMu(1), 0.0180069629759068, 1e-8);
  EXPECT_NEAR(pthermo->GetMuRatio(1), 8.18, 1e-8);
  EXPECT_NEAR(pthermo->GetInvMuRatio(1), 0.122249388753056, 1e-8);

  EXPECT_NEAR(pthermo->GetCvRatioMass(1), 0.1611002444987775, 1e-8);
  EXPECT_NEAR(pthermo->GetCvRatioMole(1), 1.3178, 1e-8);
  EXPECT_NEAR(pthermo->GetCvMassRef(1), 1521.1890586797067, 1e-8);

  EXPECT_NEAR(pthermo->GetCpRatioMass(1), 0.15, 1e-8);
  EXPECT_NEAR(pthermo->GetCpRatioMole(1), 1.227, 1e-8);
  EXPECT_NEAR(pthermo->GetCpMassRef(1), 1982.925, 1e-8);
}

TEST_F(TestThermodynamics, ammonia_vapor) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetMu(2), 0.016994346476, 1e-8);
  EXPECT_NEAR(pthermo->GetMuRatio(2), 7.72, 1e-8);
  EXPECT_NEAR(pthermo->GetInvMuRatio(2), 0.12953367875, 1e-8);

  EXPECT_NEAR(pthermo->GetCvRatioMass(2), 0.0573865284, 1e-8);
  EXPECT_NEAR(pthermo->GetCvRatioMole(2), 0.443024, 1e-8);
  EXPECT_NEAR(pthermo->GetCvMassRef(2), 541.8722953367, 1e-8);

  EXPECT_NEAR(pthermo->GetCpRatioMass(2), 0.078, 1e-8);
  EXPECT_NEAR(pthermo->GetCpRatioMole(2), 0.60216, 1e-8);
  EXPECT_NEAR(pthermo->GetCpMassRef(2), 1031.1210, 1e-8);
}

TEST_F(TestThermodynamics, water_cloud) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetCvRatioMass(3), 0.462, 1e-8);
  EXPECT_NEAR(pthermo->GetCvRatioMole(3), 3.77916, 1e-8);
  EXPECT_NEAR(pthermo->GetCvMassRef(3), 4362.435, 1e-8);

  EXPECT_NEAR(pthermo->GetCpRatioMass(3), 0.33, 1e-8);
  EXPECT_NEAR(pthermo->GetCpRatioMole(3), 2.6994, 1e-8);
  EXPECT_NEAR(pthermo->GetCpMassRef(3), 4362.435, 1e-8);

  EXPECT_NEAR(pthermo->GetLatentEnergyMass(3, 0.), 3133644.9358679, 1e-6);
  EXPECT_NEAR(pthermo->GetLatentEnergyMole(3, 0.), 56427.428339812, 1e-6);
};

TEST_F(TestThermodynamics, ammonia_cloud) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetCvRatioMass(4), 0.224, 1e-8);
  EXPECT_NEAR(pthermo->GetCvRatioMole(4), 1.72928, 1e-8);
  EXPECT_NEAR(pthermo->GetCvMassRef(4), 2115.12, 1e-8);

  EXPECT_NEAR(pthermo->GetCpRatioMass(4), 0.16, 1e-8);
  EXPECT_NEAR(pthermo->GetCpRatioMole(4), 1.2352, 1e-8);
  EXPECT_NEAR(pthermo->GetCpMassRef(4), 2115.12, 1e-8);

  EXPECT_NEAR(pthermo->GetLatentEnergyMass(4, 0.), 2262832.9904145, 1e-6);
  EXPECT_NEAR(pthermo->GetLatentEnergyMole(4, 0.), 38455.367856516, 1e-6);
};

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
