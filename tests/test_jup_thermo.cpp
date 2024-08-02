// C/C++
#include <algorithm>
#include <numeric>

// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/parameter_input.hpp>

// cantera
#include <cantera/kinetics.h>
#include <cantera/thermo.h>

// snap
#include <air_parcel.hpp>
#include <snap/thermodynamics/thermodynamics.hpp>

class TestThermodynamics : public testing::Test {
 protected:
  ParameterInput* pinput;
  char fname[80] = "/tmp/tempfile.XXXXXX";

  void CreateInputFile() {
    const char* config = R"(
<problem>
thermodynamics_config = jup-thermo.yaml
)";
    // write to file
    mkstemp(fname);
    std::ofstream outfile(fname);
    outfile << config;
  }

  virtual void SetUp() {
    CreateInputFile();

    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open(fname, IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    Thermodynamics::InitFromAthenaInput(pinput);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    delete pinput;
  }
};

TEST_F(TestThermodynamics, molecular_weight) {
  auto pthermo = Thermodynamics::GetInstance();

  EXPECT_NEAR(pthermo->GetRd(), 3571.66, 1e-2);

  Real mud = Constants::Rgas / pthermo->GetRd();

  EXPECT_NEAR(mud / pthermo->GetInvMuRatio(0), 0.0023278, 1e-3);
  EXPECT_NEAR(mud / pthermo->GetInvMuRatio(1), 0.018, 1e-3);
  EXPECT_NEAR(mud / pthermo->GetInvMuRatio(2), 0.017, 1e-3);
  EXPECT_NEAR(mud / pthermo->GetInvMuRatio(3), 0.034, 1e-3);
}

TEST_F(TestThermodynamics, vapors) {
  auto pthermo = Thermodynamics::GetInstance();

  // 1. H2O
  EXPECT_NEAR(pthermo->GetInvMuRatio(1), 0.129, 1e-3);
  EXPECT_NEAR(pthermo->GetCvRatio(1), 0.18, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatio(1), 0.166, 1e-2);

  // 2. NH3
  EXPECT_NEAR(pthermo->GetInvMuRatio(2), 0.1367, 1e-3);
  EXPECT_NEAR(pthermo->GetCvRatio(2), 0.1912, 1e-2);
  EXPECT_NEAR(pthermo->GetCpRatio(2), 0.1756, 1e-2);

  // 3. H2S
  EXPECT_NEAR(pthermo->GetInvMuRatio(3), 0.0685, 1e-3);
  EXPECT_NEAR(pthermo->GetCvRatio(3), 0.0956, 1e-3);
  EXPECT_NEAR(pthermo->GetCpRatio(3), 0.0878, 1e-3);

  // 4. H2O(l)
  EXPECT_NEAR(pthermo->GetInvMuRatio(4), 0.129, 1e-3);
  EXPECT_NEAR(pthermo->GetCvRatio(4), 0.468, 1e-3);
  EXPECT_NEAR(pthermo->GetCpRatio(4), 0.334, 1e-3);

  // 5. NH3(l)
  EXPECT_NEAR(pthermo->GetCvRatio(5), 0.526, 1e-3);
  EXPECT_NEAR(pthermo->GetCpRatio(5), 0.375, 1e-3);

  // 6. H2O(s)
  // 7. NH3(s)
  // 8. NH4SH(s)
}

TEST_F(TestThermodynamics, equilibrium_tp) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> yfrac = {0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  std::vector<Real> prim(NHYDRO, 0.);

  pthermo->SetMassFractions<Real>(yfrac.data());
  pthermo->SetTemperature(300.);
  pthermo->SetPressure(1.e5);

  pthermo->EquilibrateTP();
  pthermo->GetPrimitive<Real>(prim.data());

  // density
  EXPECT_NEAR(prim[0], 0.447323, 1e-5);
  // H2O
  EXPECT_NEAR(prim[1], 0.0571486, 1e-5);
  // NH3
  EXPECT_NEAR(prim[2], 0.333324, 1e-5);
  // H2S
  EXPECT_NEAR(prim[3], 0.166676, 1e-5);
  // H2O(l)
  EXPECT_NEAR(prim[4], 0.342851, 1e-5);
  // NH3(l)
  EXPECT_NEAR(prim[5], 0., 1e-3);
  // H2O(s)
  EXPECT_NEAR(prim[6], 0., 1e-3);
  // NH3(s)
  EXPECT_NEAR(prim[7], 0., 1e-3);
  // NH4SH(s)
  EXPECT_NEAR(prim[8], 0., 1e-3);
  // vel1
  EXPECT_NEAR(prim[9], 0., 1e-3);
  // vel2
  EXPECT_NEAR(prim[10], 0., 1e-3);
  // vel3
  EXPECT_NEAR(prim[11], 0., 1e-3);
  // pressure
  EXPECT_NEAR(prim[12], 100000, 1e-3);
}

TEST_F(TestThermodynamics, equilibrium_uv) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> yfrac = {0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  std::vector<Real> prim(NHYDRO, 0.);

  std::cout << "Before: " << std::endl;
  for (int i = 0; i < yfrac.size(); ++i) {
    std::cout << yfrac[i] << ", ";
  }
  std::cout << std::endl;

  pthermo->SetMassFractions<Real>(yfrac.data());
  pthermo->SetTemperature(300.);
  pthermo->SetPressure(1.e5);

  pthermo->EquilibrateUV();
  pthermo->GetPrimitive<Real>(prim.data());

  std::cout << "After: " << std::endl;
  for (int i = 0; i < prim.size(); ++i) {
    std::cout << prim[i] << ", ";
  }
  std::cout << std::endl;
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
