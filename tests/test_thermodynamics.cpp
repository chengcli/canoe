// C/C++
#include <algorithm>
#include <numeric>

// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// cantera
#include <cantera/kinetics.h>
#include <cantera/thermo.h>

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

class TestThermodynamics : public testing::Test {
 protected:
  ParameterInput* pinput;

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

  EXPECT_NEAR(pthermo->GetRd(), 3571.66, 1e-2);
  EXPECT_NEAR(pthermo->GetMu(0), 0.0023278, 1e-2);
}

TEST_F(TestThermodynamics, water_vapor) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetMu(1), 0.018, 1e-2);
  EXPECT_NEAR(pthermo->GetMuRatio(1), 7.739, 1e-2);
  EXPECT_NEAR(pthermo->GetInvMuRatio(1), 0.129, 1e-2);

  EXPECT_NEAR(pthermo->GetCvRatioMass(1), 0.18, 1e-2);
  EXPECT_NEAR(pthermo->GetCvRatioMole(1), 1.4, 1e-2);
  EXPECT_NEAR(pthermo->GetCvMassRef(1), 1614.5, 1e-1);

  EXPECT_NEAR(pthermo->GetCpRatioMass(1), 0.166, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatioMole(1), 1.285, 1e-2);
  EXPECT_NEAR(pthermo->GetCpMassRef(1), 2076.05, 1e-1);
}

TEST_F(TestThermodynamics, ammonia_vapor) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetMu(2), 0.017, 1e-2);
  EXPECT_NEAR(pthermo->GetMuRatio(2), 7.316, 1e-2);
  EXPECT_NEAR(pthermo->GetInvMuRatio(2), 0.1367, 1e-2);

  EXPECT_NEAR(pthermo->GetCvRatioMass(2), 0.1912, 1e-2);
  EXPECT_NEAR(pthermo->GetCvRatioMole(2), 1.4, 1e-2);
  EXPECT_NEAR(pthermo->GetCvMassRef(2), 1707.8, 1e-2);

  EXPECT_NEAR(pthermo->GetCpRatioMass(2), 0.1756, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatioMole(2), 1.285, 1e-2);
  EXPECT_NEAR(pthermo->GetCpMassRef(2), 2196, 1e-2);
}

TEST_F(TestThermodynamics, water_cloud) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetCvRatioMass(4), 0.462, 1e-2);
  EXPECT_NEAR(pthermo->GetCvRatioMole(4), 3.623, 1e-2);
  EXPECT_NEAR(pthermo->GetCvMassRef(4), 4179.85, 1e-2);

  EXPECT_NEAR(pthermo->GetCpRatioMass(4), 0.33, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatioMole(4), 2.588, 1e-2);
  EXPECT_NEAR(pthermo->GetCpMassRef(4), 4179.85, 1e-2);
};

TEST_F(TestThermodynamics, ammonia_cloud) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetCvRatioMass(5), 0.526, 1e-2);
  EXPECT_NEAR(pthermo->GetCvRatioMole(5), 3.849, 1e-2);
  EXPECT_NEAR(pthermo->GetCvMassRef(5), 4697.32, 1e-2);

  EXPECT_NEAR(pthermo->GetCpRatioMass(5), 0.375, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatioMole(5), 2.74, 1e-2);
  EXPECT_NEAR(pthermo->GetCpMassRef(5), 4697.32, 1e-2);
};

TEST_F(TestThermodynamics, equilibrium_tp) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> yfrac = {0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  auto kinetics = pthermo->Kinetics();
  auto& thermo = kinetics->thermo();

  thermo.setMassFractionsPartial(yfrac.data());
  thermo.setTemperature(300.);
  thermo.setPressure(1.e5);

  std::cout << "T = " << thermo.temperature() << std::endl;
  std::cout << "P = " << thermo.pressure() << std::endl;
  std::cout << "D = " << thermo.density() << std::endl;

  std::cout << "Number of reactions: " << kinetics->nReactions() << std::endl;

  pthermo->EquilibrateTP();
}

TEST_F(TestThermodynamics, equilibrium_uv) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> yfrac = {0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  auto kinetics = pthermo->Kinetics();
  auto& thermo = kinetics->thermo();

  thermo.setMassFractionsPartial(yfrac.data());
  thermo.setTemperature(200.);
  thermo.setPressure(1.e5);

  std::cout << "T = " << thermo.temperature() << std::endl;
  std::cout << "P = " << thermo.pressure() << std::endl;
  std::cout << "D = " << thermo.density() << std::endl;

  std::cout << "Number of reactions: " << kinetics->nReactions() << std::endl;

  pthermo->EquilibrateUV();
}

TEST_F(TestThermodynamics, saturation_adjust) {
  auto pthermo = Thermodynamics::GetInstance();

  int iH2O = 1;
  int iNH3 = 2;
  int iH2Oc = 0;
  int iNH3c = 1;

  std::vector<AirParcel> air_column(1);

  auto& air = air_column[0];
  air.SetType(AirParcel::Type::MoleFrac);

  air.SetZero();

  air.w[IDN] = 160.;
  air.w[IPR] = 7.E4;
  air.w[iH2O] = 0.02;
  air.w[iNH3] = 0.10;

  pthermo->SaturationAdjustment(air_column);

  EXPECT_NEAR(air.w[IDN], 206.41192408792, 1e-8);
  EXPECT_NEAR(air.w[IPR], 88499.534594159, 1e-8);

  Real svp1 = sat_vapor_p_H2O_BriggsS(air.w[IDN]);

  Real qgas = 1.;
#pragma omp parallel for reduction(+ : qgas)
  for (int n = 0; n < NCLOUD; ++n) qgas += -air.c[n];

  EXPECT_NEAR(air.w[iH2O] / qgas * air.w[IPR] / svp1, 1., 0.01);
  EXPECT_NEAR(air.w[iNH3], 0.1, 1e-8);

  EXPECT_NEAR(air.c[iH2Oc], 0.02 - air.w[iH2O], 1e-8);
  EXPECT_NEAR(air.c[iNH3c], 0.1 - air.w[iNH3], 1e-8);
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
