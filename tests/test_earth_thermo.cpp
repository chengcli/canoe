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
    std::vector<Real> yfrac = {0.9804, 0.0196};

    pthermo->SetMassFractions<Real>(yfrac.data());
    pthermo->EquilibrateTP(289.85, 1.e5);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    delete pinput;
  }
};

TEST_F(TestThermodynamics, equilibrate_tp) {
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<Real> prim(NHYDRO, 0.);
  pthermo->GetPrimitive<Real>(prim.data());

  for (int i = 0; i < NHYDRO; i++) {
    std::cout << prim[i] << std::endl;
  }
}

TEST_F(TestThermodynamics, cal_dlnT_dlnP) {
  auto pthermo = Thermodynamics::GetInstance();
  pthermo->Extrapolate_inplace(/*dz=*/100., /*method=*/"reversible",
                               /*grav=*/9.8);

  std::vector<Real> prim(NHYDRO, 0.);
  pthermo->GetPrimitive<Real>(prim.data());

  for (int i = 0; i < NHYDRO; i++) {
    std::cout << prim[i] << std::endl;
  }
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);  // needed for MPI initialization

  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  Application::Destroy();
  return result;
}
