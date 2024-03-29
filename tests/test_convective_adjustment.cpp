// C/C++
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// scm
#include <single_column/single_column.hpp>

// canoe
#include <air_parcel.hpp>

class TestConvectiveAdjustment : public testing::Test {
 protected:
  Mesh* pmesh;
  ParameterInput* pinput;
  Real Ps, Ts, grav;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_convective_adjustment.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    Thermodynamics::InitFromAthenaInput(pinput);

    // set up mesh
    int restart = false;
    int mesh_only = false;
    pmesh = new Mesh(pinput, mesh_only);

    // set up components
    for (int b = 0; b < pmesh->nblocal; ++b) {
      MeshBlock* pmb = pmesh->my_blocks(b);
      pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
    }

    Ps = pinput->GetReal("problem", "Ps");
    Ts = pinput->GetReal("problem", "Ts");
    grav = -pinput->GetReal("hydro", "grav_acc1");
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    Thermodynamics::Destroy();
    IndexMap::Destroy();
    delete pinput;
    delete pmesh;
  }
};

TEST_F(TestConvectiveAdjustment, Unstable) {
  auto pmb = pmesh->my_blocks(0);
  auto pthermo = Thermodynamics::GetInstance();
  auto phydro = pmb->phydro;
  auto pcoord = pmb->pcoord;

  int ks = pmb->ks;
  int js = pmb->js;
  int is = pmb->is;

  AirColumn ac(8);

  ac[0].ToMoleFraction();
  ac[0].w[IPR] = 9.61664e+06;
  ac[0].w[IDN] = 4929.37;
  // theta = 4984.74

  ac[1].ToMoleFraction();
  ac[1].w[IPR] = 8.87938e+06;
  ac[1].w[IDN] = 4819.53;
  // theta = 4986

  ac[2].ToMoleFraction();
  ac[2].w[IPR] = 8.18359e+06;
  ac[2].w[IDN] = 4717.91;
  // theta = 4996

  ac[3].ToMoleFraction();
  ac[3].w[IPR] = 7.52891e+06;
  ac[3].w[IDN] = 4592.92;
  // theta = 4980.9

  ac[4].ToMoleFraction();
  ac[4].w[IPR] = 6.91071e+06;
  ac[4].w[IDN] = 4467.53;
  // theta = 4964.98

  ac[5].ToMoleFraction();
  ac[5].w[IPR] = 6.32785e+06;
  ac[5].w[IDN] = 4354.14;
  // theta = 4962.34

  ac[6].ToMoleFraction();
  ac[6].w[IPR] = 5.78069e+06;
  ac[6].w[IDN] = 4237.51;
  // theta = 4955.84

  ac[7].ToMoleFraction();
  ac[7].w[IPR] = 5.26754e+06;
  ac[7].w[IDN] = 4139.14;
  // theta = 4971.08

  Real mass0 = 0., enthalpy0 = 0.;
  for (int i = 0; i <= 7; ++i) {
    auto air = ac[i].ToMassConcentration();
    Real density = ac[i].w[IDN];
    mass0 += density;
    enthalpy0 += (ac[i].w[IEN] + density * grav * pcoord->x1v(i));
  }

  pmb->pimpl->pscm->ConvectiveAdjustment(ac, pmb->ks, pmb->js, 0, 7);

  Real mass1 = 0., enthalpy1 = 0.;
  for (int i = 0; i <= 7; ++i) {
    auto air = ac[i];
    air.ToMassConcentration();

    Real density = ac[i].w[IDN];
    mass1 += density;
    enthalpy1 += (ac[i].w[IEN] + density * grav * pcoord->x1v(i));
  }

  for (int i = 0; i <= 7; ++i) {
    auto air = ac[i];
    air.ToMoleFraction();
    std::cout << i << "," << air.w[IPR] << "," << air.w[IDN] << ","
              << air.theta(Ps) << std::endl;
  }

  EXPECT_NEAR(mass1 / mass0, 1., 1.e-2);
  EXPECT_NEAR(enthalpy1 / enthalpy0, 1., 1.e-2);
}

TEST_F(TestConvectiveAdjustment, NegativeEnergy) {
  auto pmb = pmesh->my_blocks(0);
  auto pthermo = Thermodynamics::GetInstance();
  auto phydro = pmb->phydro;
  auto pcoord = pmb->pcoord;

  int ks = pmb->ks;
  int js = pmb->js;
  int is = pmb->is;

  AirColumn ac(8);

  ac[0].ToMoleFraction();
  ac[0].w[IPR] = 9.61664e+06;
  ac[0].w[IDN] = 4929.37;
  // theta = 4984.74

  ac[1].ToMoleFraction();
  ac[1].w[IPR] = 8.87938e+06;
  ac[1].w[IDN] = 4819.53;
  // theta = 4986

  ac[2].ToMoleFraction();
  ac[2].w[IPR] = -1.e4;
  ac[2].w[IDN] = 4717.91;
  // theta = Nan

  ac[3].ToMoleFraction();
  ac[3].w[IPR] = 7.52891e+06;
  ac[3].w[IDN] = 4592.92;
  // theta = 4980.9

  ac[4].ToMoleFraction();
  ac[4].w[IPR] = 6.91071e+06;
  ac[4].w[IDN] = 4467.53;
  // theta = 4964.98

  Real mass0 = 0., enthalpy0 = 0.;
  for (int i = 0; i <= 4; ++i) {
    auto air = ac[i];
    air.ToMassConcentration();

    Real density = air.w[IDN];
    mass0 += density;
    enthalpy0 += (air.w[IEN] + density * grav * pcoord->x1v(i));

    air.ToMoleFraction();
    std::cout << i << "," << air.w[IPR] << "," << air.w[IDN] << ","
              << air.theta(Ps) << std::endl;
  }

  pmb->pimpl->pscm->ConvectiveAdjustment(ac, pmb->ks, pmb->js, 0, 4);

  Real mass1 = 0., enthalpy1 = 0.;
  for (int i = 0; i <= 4; ++i) {
    auto air = ac[i];
    air.ToMassConcentration();

    Real density = air.w[IDN];
    mass1 += density;
    enthalpy1 += (air.w[IEN] + density * grav * pcoord->x1v(i));

    air.ToMoleFraction();
    std::cout << i << "," << air.w[IPR] << "," << air.w[IDN] << ","
              << air.theta(Ps) << std::endl;
  }

  EXPECT_NEAR(mass1 / mass0, 1., 1.e-3);
  EXPECT_NEAR(enthalpy1 / enthalpy0, 1., 1.e-3);
}

TEST_F(TestConvectiveAdjustment, RandomProfile) {
  return;

  auto pmb = pmesh->my_blocks(0);
  auto pcoord = pmb->pcoord;

  // prepare random real number
  std::random_device rd;
  std::mt19937 gen(rd());  // Mersenne Twister generator

  // Define the distribution to be uniform for real numbers
  // between 0.0 and 1.0
  std::uniform_real_distribution<> dis(-20., 20.0);

  // output file
  std::ofstream outFile1("ac_before.csv");
  outFile1 << "i,T_fluct,x1,pres,temp,theta" << std::endl;

  // construt the air column
  std::vector<AirParcel> vector_ac;

  AirParcel air(AirParcel::Type::MoleFrac);
  air.SetZero();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;

  auto pthermo = Thermodynamics::GetInstance();

  // half a grid to cell center
  pthermo->Extrapolate(&air, pcoord->dx1f(pmb->is) / 2., "reversible", grav);

  AirColumn ac(pmb->ncells1);

  int js = pmb->js, ks = pmb->ks;
  for (int i = pmb->is; i <= pmb->ie; ++i) {
    Real T_fluct = dis(gen);
    air.ToMoleFraction();
    air.w[IDN] += T_fluct;

    outFile1 << i << "," << T_fluct << "," << pcoord->x1v(i) << ","
             << air.w[IPR] << "," << air.w[IDN] << "," << air.theta(Ps)
             << std::endl;

    ac[i] = air;
    pthermo->Extrapolate(&air, pcoord->dx1f(i), "pseudo", grav);
  }
  outFile1.close();

  auto pscm = pmb->pimpl->pscm;
  std::vector<std::array<int, 2>> ranges;
  pscm->FindUnstableRange(ac, pmb->is, pmb->ie, ranges);

  for (auto range : ranges) {
    std::cout << range[0] << " " << range[1] << std::endl;
    pscm->ConvectiveAdjustment(ac, pmb->ks, pmb->js, range[0], range[1]);
  }

  // output air column after convective adjustment
  std::ofstream outFile2("ac_after.csv");
  outFile2 << "i,x1,pres,temp,theta" << std::endl;

  for (int i = pmb->is; i <= pmb->ie; ++i) {
    air = ac[i];
    air.ToMoleFraction();
    outFile2 << i << "," << pcoord->x1v(i) << "," << air.w[IPR] << ","
             << air.w[IDN] << "," << air.theta(Ps) << std::endl;
  }
  outFile2.close();
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
