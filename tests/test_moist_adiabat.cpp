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

// microphysics
#include <microphysics/microphysics.hpp>

double sat_vapor_p_H2O(double T) {
  double betal = 24.845, gammal = 4.986009, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

// water svp
void enroll_vapor_function_H2O(Thermodynamics::SVPFunc1Container& svp_func1) {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1[iH2O][0] = [](AirParcel const& qfrac, int, int) {
    return sat_vapor_p_H2O(qfrac.w[IDN]);
  };
}

void Thermodynamics::enrollVaporFunctions() {
  enroll_vapor_function_H2O(svp_func1_);
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
  Mesh* pmesh;
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
    delete pmesh;
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

  AirParcel air(AirParcel::Type::MassFrac);
  air.w[iH2O] = qt;
  air.w[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.w[iH2Oc] = 0.;

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

  AirParcel air(AirParcel::Type::MassFrac);
  air.w[iH2O] = qt;
  air.w[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.w[iH2Oc] = 0.;

  pthermo->Extrapolate(&air, dz, "reversible", grav);

  EXPECT_NEAR(air.w[IDN], 289.392118923, 1e-8);
  EXPECT_NEAR(air.w[IPR], 98825.8592854826, 1e-8);
  EXPECT_NEAR(air.w[iH2O], 0.0184640437464, 1e-8);
  EXPECT_NEAR(air.w[iH2Oc], 0.0127248713312, 1e-8);
}

TEST_F(TestMoistAdiabat, moist_static_energy) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmesh->my_blocks(0);

  AirParcel air(AirParcel::Type::MassFrac);
  air.w[iH2O] = qt;
  air.w[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.w[iH2Oc] = 0.;

  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);
  air.w[iH2O] += rates[0];
  air.w[iH2Oc] += rates[1];

  // first grid
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is, air);

  Real mse1 = pthermo->MoistStaticEnergy(pmb, 0., ks, js, is);
  EXPECT_NEAR(mse1, 272872.16946, 1.E-4);

  // second grid
  pthermo->Extrapolate(&air, dz, "reversible", grav);
  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is + 1, air);

  Real mse2 = pthermo->MoistStaticEnergy(pmb, grav * dz, ks, js, is + 1);
  EXPECT_NEAR(mse2, 272872.16971, 1.E-4);
}

TEST_F(TestMoistAdiabat, equivalent_potential_temp) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmesh->my_blocks(0);

  AirParcel air(AirParcel::Type::MassFrac);
  air.w[iH2O] = qt;
  air.w[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.w[iH2Oc] = 0.;

  auto rates = pthermo->TryEquilibriumTP_VaporCloud(air, iH2O);
  air.w[iH2O] += rates[0];
  air.w[iH2Oc] += rates[1];

  // first grid
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is, air);

  Real theta_e1 = pthermo->EquivalentPotentialTemp(pmb, Ps, iH2O, ks, js, is);
  EXPECT_NEAR(theta_e1, 320.54669065133, 1.E-8);

  // second grid
  pthermo->Extrapolate(&air, dz, "reversible", grav);
  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is + 1, air);

  Real theta_e2 =
      pthermo->EquivalentPotentialTemp(pmb, Ps, iH2O, ks, js, is + 1);
  EXPECT_NEAR(theta_e2, 320.5466916420, 1.E-4);
}

TEST_F(TestMoistAdiabat, saturation_adjustment) {
  auto pthermo = Thermodynamics::GetInstance();
  auto pmb = pmesh->my_blocks(0);
  Real vx = 100.;
  Real vy = 200.;
  Real vz = 300.;

  AirColumn air_column(1);

  auto& air = air_column[0];
  air.SetType(AirParcel::Type::MassFrac);

  air.w[iH2O] = qt;
  air.w[iH2Oc] = 0.;

  air.ToMoleFraction();
  air.w[IPR] = Ps;
  air.w[IDN] = Ts;
  air.w[IVX] = vx;
  air.w[IVY] = vy;
  air.w[IVZ] = vz;
  air.w[iH2Oc] = 0.;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  pthermo->SaturationAdjustment(air_column);

  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is, air);
  AirParcelHelper::distribute_to_conserved(pmb, ks, js, is, air);

  air.ToMassFraction();
  EXPECT_NEAR(air.w[IDN], 1.187901949988, 1e-8);
  EXPECT_NEAR(air.w[IPR], 101934.372666, 1e-6);
  EXPECT_NEAR(air.w[iH2O] + air.w[iH2Oc], 0.0196, 1e-8);
  air.ToMoleFraction();

  Real drho = pmb->pimpl->pmicro->u(0, ks, js, is) / 2.;
  Real cv = pthermo->GetCvMassRef(iH2O);
  Real temp = air.w[IDN];
  pmb->phydro->u(iH2O, ks, js, is) += drho;

  pmb->phydro->u(IVX, ks, js, is) += drho * vx;
  pmb->phydro->u(IVY, ks, js, is) += drho * vy;
  pmb->phydro->u(IVZ, ks, js, is) += drho * vz;
  pmb->phydro->u(IEN, ks, js, is) +=
      0.5 * drho * (vx * vx + vy * vy + vz * vz) + drho * cv * temp;

  pmb->pimpl->pmicro->u(0, ks, js, is) -= drho;
  air = AirParcelHelper::gather_from_conserved(pmb, ks, js, is);
  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is, air);

  Real theta_e1 = pthermo->EquivalentPotentialTemp(pmb, Ps, iH2O, ks, js, is);

  pthermo->SaturationAdjustment(air_column);

  AirParcelHelper::distribute_to_primitive(pmb, ks, js, is, air);
  Real theta_e2 = pthermo->EquivalentPotentialTemp(pmb, Ps, iH2O, ks, js, is);

  EXPECT_NEAR(theta_e1, 343.73446046, 1e-8);
  EXPECT_NEAR(theta_e2, 343.73573051, 1e-8);
}

int main(int argc, char* argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
