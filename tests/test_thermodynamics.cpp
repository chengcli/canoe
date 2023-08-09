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

class TestThermodynamics : public testing::Test {
 protected:
  ParameterInput *pinput;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_thermodynamics.inp", IOWrapper::FileMode::read);

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

TEST_F(TestThermodynamics, latent_heats) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> rates(3);
  rates[0] = -1.;
  rates[1] = 0.5;
  rates[2] = 0.5;

  EXPECT_NEAR(pthermo->GetLatentEnergyMole(3, 300.), 43573.103798572, 1e-6);
  EXPECT_NEAR(pthermo->GetLatentHeatMole(1, rates, 300.), 46067.442398572,
              1e-6);

  EXPECT_NEAR(pthermo->GetLatentEnergyMass(3, 300.), 2419791.93586797, 1e-6);
  EXPECT_NEAR(pthermo->GetLatentHeatMass(1, rates, 300.), 2419791.9358679,
              1e-6);
}

TEST_F(TestThermodynamics, equilibrium_water) {
  auto pthermo = Thermodynamics::GetInstance();

  int iH2O = 1;
  int iNH3 = 2;

  Variable air(Variable::Type::MoleFrac);
  air.SetZero();

  air.w[IDN] = 300.;
  air.w[IPR] = 7.E5;
  air.w[iH2O] = 0.2;
  air.w[iNH3] = 0.1;

  // water
  Real svp = sat_vapor_p_H2O_BriggsS(air.w[IDN]);
  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);

  EXPECT_NEAR(rates[0], svp / air.w[IPR] - air.w[iH2O], 1e-3);
  EXPECT_NEAR(rates[1], 0.19592911846053, 1e-8);
  EXPECT_NEAR(rates[2], 0.0, 1e-8);
}

TEST_F(TestThermodynamics, equilibrium_ammonia) {
  auto pthermo = Thermodynamics::GetInstance();

  int iH2O = 1;
  int iNH3 = 2;

  Variable air(Variable::Type::MoleFrac);
  air.SetZero();

  air.w[IDN] = 160.;
  air.w[IPR] = 7.E4;
  air.w[iH2O] = 0.02;
  air.w[iNH3] = 0.01;

  // ammonia
  Real svp = sat_vapor_p_NH3_BriggsS(air.w[IDN]);
  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iNH3);

  EXPECT_NEAR(rates[0], svp / air.w[IPR] - air.w[iNH3], 1e-3);
  EXPECT_NEAR(rates[1], 0.0088517945331865, 1e-8);
  EXPECT_NEAR(rates[2], 0.0, 1e-8);
}

TEST_F(TestThermodynamics, saturation_adjust) {
  auto pthermo = Thermodynamics::GetInstance();

  int iH2O = 1;
  int iNH3 = 2;
  int iH2Oc = 0;
  int iNH3c = 1;

  std::vector<Variable> air_column(1);

  auto &air = air_column[0];
  air.SetType(Variable::Type::MoleFrac);

  air.SetZero();

  air.w[IDN] = 160.;
  air.w[IPR] = 7.E4;
  air.w[iH2O] = 0.02;
  air.w[iNH3] = 0.10;

  pthermo->SaturationAdjustment(air_column);

  EXPECT_NEAR(air.w[IDN], 206.41191698346, 1e-8);
  EXPECT_NEAR(air.w[IPR], 88499.531838818, 1e-8);

  Real svp1 = sat_vapor_p_H2O_BriggsS(air.w[IDN]);

  Real qgas = 1.;
#pragma omp parallel for reduction(+ : qgas)
  for (int n = 0; n < NCLOUD; ++n) qgas += -air.c[n];

  EXPECT_NEAR(air.w[iH2O] / qgas * air.w[IPR] / svp1, 1., 0.01);
  EXPECT_NEAR(air.w[iNH3], 0.1, 1e-8);

  EXPECT_NEAR(air.c[iH2Oc], 0.02 - air.w[iH2O], 1e-8);
  EXPECT_NEAR(air.c[iNH3c], 0.1 - air.w[iNH3], 1e-8);
}

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
