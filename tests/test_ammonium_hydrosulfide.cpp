// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// canoe
#include <impl.hpp>
#include <index_map.hpp>
#include <variable.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/vapors/ammonia_vapors.hpp>
#include <snap/thermodynamics/vapors/water_vapors.hpp>

class TestAmmoniumHydrosulfide : public testing::Test {
 protected:
  ParameterInput *pinput;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_ammonium_hydrosulfide.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    IndexMap::InitFromAthenaInput(pinput);
    Thermodynamics::InitFromAthenaInput(pinput);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    IndexMap::Destroy();
    delete pinput;
  }
};

TEST_F(TestAmmoniumHydrosulfide, molecule) {
  auto pindex = IndexMap::GetInstance();
  auto pthermo = Thermodynamics::GetInstance();

  int iNH4SH = pindex->GetCloudId("NH4SH(s)");
  int j = iNH4SH + 1 + NVAPOR;

  EXPECT_DOUBLE_EQ(pthermo->GetMuRatio(j), 0.051111 / pthermo->GetMu(0));
  EXPECT_DOUBLE_EQ(pthermo->GetMu(j), 0.051111);
  EXPECT_DOUBLE_EQ(pthermo->GetInvMu(j), 1. / 0.051111);

  EXPECT_DOUBLE_EQ(pthermo->GetCpMassRef(j), 0.);
  EXPECT_DOUBLE_EQ(pthermo->GetCpRatioMass(j), 0.);
  EXPECT_DOUBLE_EQ(pthermo->GetCpRatioMole(j), 0.);

  EXPECT_DOUBLE_EQ(pthermo->GetCvMassRef(j), 0.);
  EXPECT_DOUBLE_EQ(pthermo->GetCvRatioMass(j), 0.);
  EXPECT_DOUBLE_EQ(pthermo->GetCvRatioMole(j), 0.);
}

TEST_F(TestAmmoniumHydrosulfide, equilibrium) {}

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
