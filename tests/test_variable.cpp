// C/C++
#include <algorithm>
#include <cmath>

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

class TestVariable : public testing::Test {
 protected:
  ParameterInput *pinput;
  Variable a;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_thermodynamics.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    IndexMap::InitFromAthenaInput(pinput);
    Thermodynamics::InitFromAthenaInput(pinput);

    a.SetType(Variable::Type::MoleFrac);

    std::fill(a.w, a.w + Variable::Size,
              1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));
    std::fill(a.x, a.x + NTRACER, 2.);

    a.w[IDN] = 1.;
    a.w[IPR] = 10.;
    a.w[IVX] = 0.1;
    a.w[IVY] = 0.2;
    a.w[IVZ] = 0.3;
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
    //
    Thermodynamics::Destroy();
    IndexMap::Destroy();
    delete pinput;
  }
};

TEST(variable, copy_constructer) {
  Variable a(Variable::Type::MoleFrac);

  std::fill(a.w, a.w + Variable::Size, 1.0);

  Variable b(a);
  std::fill(b.w, b.w + Variable::Size, 2.0);

  EXPECT_DOUBLE_EQ(a.w[0], 1.0);
  EXPECT_DOUBLE_EQ(a.w[1], 1.0);
  EXPECT_DOUBLE_EQ(a.w[2], 1.0);

  EXPECT_DOUBLE_EQ(b.w[0], 2.0);
  EXPECT_DOUBLE_EQ(b.w[1], 2.0);
  EXPECT_DOUBLE_EQ(b.w[2], 2.0);
}

TEST(variable, assignment_operator) {
  Variable a(Variable::Type::MoleFrac);

  std::fill(a.w, a.w + Variable::Size, 1.0);

  Variable b(Variable::Type::MoleFrac);
  b = a;
  std::fill(b.w, b.w + Variable::Size, 2.0);

  EXPECT_DOUBLE_EQ(a.w[0], 1.0);
  EXPECT_DOUBLE_EQ(a.w[1], 1.0);
  EXPECT_DOUBLE_EQ(a.w[2], 1.0);

  EXPECT_DOUBLE_EQ(b.w[0], 2.0);
  EXPECT_DOUBLE_EQ(b.w[1], 2.0);
  EXPECT_DOUBLE_EQ(b.w[2], 2.0);
}

TEST_F(TestVariable, MoleFraction_MassFraction) {
  // mole fraction to mass fraction
  a.ToMassFraction();
  a.ToMoleFraction();

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_DOUBLE_EQ(a.w[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_DOUBLE_EQ(a.c[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  for (int n = 0; n < NCHEMISTRY; ++n)
    EXPECT_DOUBLE_EQ(a.q[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  EXPECT_NEAR(a.w[IDN], 1., 1e-10);
  EXPECT_DOUBLE_EQ(a.w[IPR], 10.);
  EXPECT_DOUBLE_EQ(a.w[IVX], 0.1);
  EXPECT_DOUBLE_EQ(a.w[IVY], 0.2);
  EXPECT_DOUBLE_EQ(a.w[IVZ], 0.3);
}

TEST_F(TestVariable, MoleFraction_MassConcentration) {
  // mole fraction to mass concentration
  a.ToMassConcentration();
  a.ToMoleFraction();

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_NEAR(a.w[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_NEAR(a.c[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  for (int n = 0; n < NCHEMISTRY; ++n)
    EXPECT_NEAR(a.q[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  EXPECT_DOUBLE_EQ(a.w[IDN], 1.);
  EXPECT_NEAR(a.w[IPR], 10., 1e-10);
  EXPECT_DOUBLE_EQ(a.w[IVX], 0.1);
  EXPECT_DOUBLE_EQ(a.w[IVY], 0.2);
  EXPECT_DOUBLE_EQ(a.w[IVZ], 0.3);
}

TEST_F(TestVariable, MassFraction_MassConcentration) {
  a.SetType(Variable::Type::MassFrac);

  // mass fraction to mass concentration
  a.ToMassConcentration();
  a.ToMassFraction();

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_DOUBLE_EQ(a.w[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_DOUBLE_EQ(a.c[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  for (int n = 0; n < NCHEMISTRY; ++n)
    EXPECT_DOUBLE_EQ(a.q[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY));

  EXPECT_DOUBLE_EQ(a.w[IDN], 1.);
  EXPECT_DOUBLE_EQ(a.w[IPR], 10.);
  EXPECT_DOUBLE_EQ(a.w[IVX], 0.1);
  EXPECT_DOUBLE_EQ(a.w[IVY], 0.2);
  EXPECT_DOUBLE_EQ(a.w[IVZ], 0.3);
}

TEST_F(TestVariable, MoleFraction_MoleConcentration) {
  a.ToMoleConcentration();
  a.ToMoleFraction();

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_NEAR(a.w[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_NEAR(a.c[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  for (int n = 0; n < NCHEMISTRY; ++n)
    EXPECT_NEAR(a.q[n], 1. / (1 + NVAPOR + NCLOUD + NCHEMISTRY), 1.e-10);

  EXPECT_NEAR(a.w[IDN], 1., 1e-10);
  EXPECT_NEAR(a.w[IPR], 10., 1e-10);
  EXPECT_DOUBLE_EQ(a.w[IVX], 0.1);
  EXPECT_DOUBLE_EQ(a.w[IVY], 0.2);
  EXPECT_DOUBLE_EQ(a.w[IVZ], 0.3);
}

TEST_F(TestVariable, MassFraction_MoleConcentration) {}

TEST_F(TestVariable, MassConcentration_MoleConcentration) {}

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
