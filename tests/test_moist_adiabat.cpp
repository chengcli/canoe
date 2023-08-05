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

// water svp
void Thermodynamics::enrollVaporFunctionH2O() {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1_[iH2O][0] = [](Variable const& qfrac, int, int) {
    return sat_vapor_p_H2O_liquid_Ideal(qfrac.w[IDN]);
  };
}

// expose private members for testing
class ThermodynamicsTestOnly : public Thermodynamics {
 public:
  std::array<Real, Size> const& GetBeta() const { return beta_; }
  std::array<Real, Size> const& GetDelta() const { return delta_; }
  std::array<Real, 1 + NVAPOR> const& GetT3() const { return t3_; }
  std::array<Real, 1 + NVAPOR> const& GetP3() const { return p3_; }
};

class TestMoistAdiabat : public testing::Test {
 protected:
  ParameterInput* pinput;
  Real Ps, Ts, qt, grav, dz;
  int iH2O, iH2Oc;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_moist_adiabat.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    auto pindex = IndexMap::InitFromAthenaInput(pinput);
    Thermodynamics::InitFromAthenaInput(pinput);

    Ps = pinput->GetReal("problem", "Ps");
    Ts = pinput->GetReal("problem", "Ts");
    qt = pinput->GetReal("problem", "qt");
    grav = -pinput->GetReal("hydro", "grav_acc1");
    dz = 100.;  // grid spacing 100 m

    iH2O = pindex->GetSpeciesId("vapor.H2O");
    iH2Oc = pindex->GetSpeciesId("cloud.H2O(c)");
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
  auto& t3 = pthermo_test->GetT3();
  auto& p3 = pthermo_test->GetP3();

  EXPECT_NEAR(beta[0], 0., 1e-8);
  EXPECT_NEAR(beta[1], 0., 1e-8);
  EXPECT_NEAR(beta[2], 24.845, 1e-8);

  EXPECT_NEAR(delta[0], 0., 1e-8);
  EXPECT_NEAR(delta[1], 0., 1e-8);
  EXPECT_NEAR(delta[2], 4.986009, 1e-8);

  EXPECT_DOUBLE_EQ(t3[0], 0.);
  EXPECT_DOUBLE_EQ(t3[1], 273.16);
  EXPECT_DOUBLE_EQ(p3[0], 0.);
  EXPECT_DOUBLE_EQ(p3[1], 611.7);

  EXPECT_NEAR(pthermo->GetLatentEnergyMass(0), 0., 1e-8);
  EXPECT_NEAR(pthermo->GetLatentEnergyMass(1), 0., 1e-8);
  EXPECT_NEAR(pthermo->GetLatentEnergyMass(2), 3136508.0151368757, 1e-8);
}

TEST_F(TestMoistAdiabat, concentration) {
  auto pthermo = Thermodynamics::GetInstance();

  Variable air(Variable::Type::MassFrac);
  air.w[iH2O] = qt;
  air.c[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.c[iH2Oc] = 0.;

  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);
  air.w[iH2O] += rates[0];
  air.w[iH2Oc] += rates[1];

  EXPECT_NEAR(air.w[IDN], Ts, 1e-8);
  EXPECT_NEAR(air.w[IPR], Ps, 1e-8);
  EXPECT_NEAR(air.w[iH2O], 0.0187935176958, 1e-8);
  EXPECT_NEAR(air.w[iH2Oc], 0.0123953973819, 1e-8);

  air.ToMassFraction();

  EXPECT_NEAR(air.w[IDN], 1.2028112737, 1e-8);
  EXPECT_NEAR(air.w[IPR], Ps, 1e-8);
  EXPECT_NEAR(air.w[iH2O], 0.0118103802559, 1e-8);
  EXPECT_NEAR(air.w[iH2Oc], 0.007789619744086, 1e-8);
}

TEST_F(TestMoistAdiabat, moist_adiabat) {
  auto pthermo = Thermodynamics::GetInstance();

  Variable air(Variable::Type::MassFrac);
  air.w[iH2O] = qt;
  air.c[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.c[iH2Oc] = 0.;

  pthermo->Extrapolate(&air, dz, Thermodynamics::Method::ReversibleAdiabat,
                       grav);

  EXPECT_NEAR(air.w[IDN], 289.392118923, 1e-8);
  EXPECT_NEAR(air.w[IPR], 98825.8592854826, 1e-8);
  EXPECT_NEAR(air.w[iH2O], 0.0184640437464, 1e-8);
  EXPECT_NEAR(air.w[iH2Oc], 0.0127248713312, 1e-8);
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
