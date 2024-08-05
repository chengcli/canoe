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
#include <snap/thermodynamics/thermodynamics.hpp>

class TestThermodynamics : public testing::Test {
 protected:
  ParameterInput* pinput;
  char fname[80] = "/tmp/tempfile.XXXXXX";

  void CreateInputFile() {
    const char* config = R"(
<problem>
thermodynamics_config = earth-thermo.yaml
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

    auto pthermo = Thermodynamics::InitFromAthenaInput(pinput);
    std::vector<Real> yfrac = {0.9804, 0.0196, 0.};

    pthermo->SetMassFractions<Real>(yfrac.data());
    pthermo->EquilibrateTP(273.16, 1.e5);

    /*auto kinetics = get_kinetics_object(pthermo);
    auto& thermo = kinetics->thermo();
    for (size_t i = 0.; i < Thermodynamics::Size; ++i) {
      auto& tp = thermo.species(i)->thermo;
      size_t n;
      int type;
      Real tlow, thigh, pref;
      std::vector<Real> coeffs(tp->nCoeffs());
      tp->reportParameters(n, type, tlow, thigh, pref, coeffs.data());
      std::cout << "type = " << type << std::endl;
      std::cout << "tlow = " << tlow << std::endl;
      std::cout << "thigh = " << thigh << std::endl;
      std::cout << "pref = " << pref << std::endl;
      std::cout << "t0 = " << coeffs[0] << std::endl;
      std::cout << "h0 = " << coeffs[1] << std::endl;
      std::cout << "s0 = " << coeffs[2] << std::endl;
      std::cout << "cp0 = " << coeffs[3] << std::endl;
    }*/
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    delete pinput;
  }
};

TEST_F(TestThermodynamics, thermo) {
  auto kinetics = get_kinetics_object(Thermodynamics::GetInstance());
  auto& thermo = kinetics->thermo();

  std::array<Real, Thermodynamics::Size> cp;
  std::array<Real, Thermodynamics::Size> enthalpy;

  thermo.getPartialMolarCp(cp.data());
  thermo.getEnthalpy_RT(enthalpy.data());

  for (int i = 0; i < Thermodynamics::Size; i++) {
    cp[i] /= thermo.molecularWeight(i);
  }

  // 1. H2O
  EXPECT_NEAR(thermo.molecularWeight(1) / thermo.molecularWeight(0), 0.621,
              1e-3);
  EXPECT_NEAR(cp[1] / cp[0], 1.166, 1e-3);

  // 2. H2O(l)
  EXPECT_NEAR(thermo.molecularWeight(2) / thermo.molecularWeight(0), 0.621,
              1e-3);
  EXPECT_NEAR(cp[2] / cp[0], 3.457, 1e-3);

  // 3. beta
  EXPECT_NEAR(enthalpy[1] - enthalpy[2], 19.858, 1e-3);
}

TEST_F(TestThermodynamics, equilibrate_tp) {
  auto pthermo = Thermodynamics::GetInstance();

  pthermo->EquilibrateTP(289.85, 1.e5);

  std::vector<Real> prim(NHYDRO, 0.);
  pthermo->GetPrimitive<Real>(prim.data());

  EXPECT_NEAR(prim[IDN], 1.2028, 1e-4);
  EXPECT_NEAR(prim[1], 0.01183, 1e-4);
  EXPECT_NEAR(prim[2], 0.00777, 1e-4);
  EXPECT_NEAR(prim[IPR], 1.e5, 1e-8);
}

TEST_F(TestThermodynamics, moist_adiabat) {
  auto pthermo = Thermodynamics::GetInstance();

  pthermo->EquilibrateTP(289.85, 1.e5);
  pthermo->Extrapolate_inplace(/*dz=*/100., /*method=*/"reversible",
                               /*grav=*/9.81);

  auto kinetics = get_kinetics_object(pthermo);
  auto& thermo = kinetics->thermo();

  EXPECT_NEAR(thermo.temperature(), 289.392, 1e-3);
  EXPECT_NEAR(thermo.pressure(), 98825.8, 0.1);
  EXPECT_NEAR(thermo.moleFraction(1), 0.01846, 1e-4);
  EXPECT_NEAR(thermo.moleFraction(2), 0.01268, 1e-4);
}

TEST_F(TestThermodynamics, equilibrate_uv) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> yfrac = {0.9804, 0.0196, 0.};
  pthermo->SetMassFractions<Real>(yfrac.data());

  auto kinetics = get_kinetics_object(pthermo);
  auto& thermo = kinetics->thermo();
  thermo.setState_TPY(289.85, 1.e5, yfrac.data());

  pthermo->EquilibrateUV();

  std::vector<Real> prim(NHYDRO, 0.);
  pthermo->GetPrimitive<Real>(prim.data());
  EXPECT_NEAR(prim[IDN], 1.187901949988, 1e-4);
  EXPECT_NEAR(prim[IPR], 101928., 0.1);
  EXPECT_NEAR(prim[1] + prim[2], 0.0196, 1e-8);
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);  // needed for MPI initialization

  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  Application::Destroy();
  return result;
}
