// C/C++
#include <algorithm>

// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>
#include <snap/thermodynamics/vapors/ammonia_vapors.hpp>
#include <snap/thermodynamics/vapors/water_vapors.hpp>

// special includes
#include <special/giants_enroll_vapor_functions_v1.hpp>

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

TEST_F(TestAmmoniumHydrosulfide, equilibrium) {
  auto pindex = IndexMap::GetInstance();
  auto pthermo = Thermodynamics::GetInstance();

  int iNH3 = pindex->GetVaporId("NH3");
  int iH2S = pindex->GetVaporId("H2S");

  AirParcel qfrac(AirParcel::Type::MoleFrac);
  qfrac.SetZero();

  qfrac.w[IDN] = 160.;
  qfrac.w[IPR] = 7.E4;
  qfrac.w[iNH3] = 0.02;
  qfrac.w[iH2S] = 0.01;

  IndexPair ij(std::minmax(iNH3, iH2S));

  auto rates = pthermo->TryEquilibriumTP_VaporVaporCloud(qfrac, ij);

  EXPECT_NEAR(rates[0], -0.01, 1e-6);
  EXPECT_NEAR(rates[1], -0.01, 1e-6);
  EXPECT_NEAR(rates[2], 0.01, 1e-6);
}

TEST_F(TestAmmoniumHydrosulfide, equilibrium2) {
  auto pindex = IndexMap::GetInstance();
  auto pthermo = Thermodynamics::GetInstance();

  int iNH3 = pindex->GetVaporId("NH3");
  int iH2S = pindex->GetVaporId("H2S");

  AirParcel qfrac(AirParcel::Type::MoleFrac);
  qfrac.SetZero();

  qfrac.w[IDN] = 250.;
  qfrac.w[IPR] = 1.E5;
  qfrac.w[iNH3] = 0.01;
  qfrac.w[iH2S] = 0.02;

  IndexPair ij(std::minmax(iNH3, iH2S));

  auto rates = pthermo->TryEquilibriumTP_VaporVaporCloud(qfrac, ij);

  EXPECT_NEAR(rates[0], -0.00376944377451, 1e-6);
  EXPECT_NEAR(rates[1], -0.00376944377451, 1e-6);
  EXPECT_NEAR(rates[2], 0.00376944377451, 1e-6);
}

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
