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
#include <snap/thermodynamics/vapors/water_vapors.hpp>

// expose private members for testing
class ThermodynamicsTestOnly : public Thermodynamics {
 public:
  std::array<Real, Size> const& GetBeta() const { return beta_; }
  std::array<Real, Size> const& GetDelta() const { return delta_; }
};

class TestMoistAdiabat : public testing::Test {
 protected:
  ParameterInput* pinput;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_moist_adiabat.inp", IOWrapper::FileMode::read);

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

TEST_F(TestMoistAdiabat, parameter) {
  auto pthermo = Thermodynamics::GetInstance();

  ThermodynamicsTestOnly const* pthermo_test =
      static_cast<ThermodynamicsTestOnly const*>(pthermo);

  auto& beta = pthermo_test->GetBeta();
  auto& delta = pthermo_test->GetDelta();

  EXPECT_NEAR(beta[0], 0., 1e-8);
  EXPECT_NEAR(beta[1], 0., 1e-8);
  EXPECT_NEAR(beta[2], 24.845, 1e-8);

  EXPECT_NEAR(delta[0], 0., 1e-8);
  EXPECT_NEAR(delta[1], 0., 1e-8);
  EXPECT_NEAR(delta[2], 4.986009, 1e-8);
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
