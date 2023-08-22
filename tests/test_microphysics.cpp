// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <impl.hpp>
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

class TestMicrophysics : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;

  virtual void SetUp() {
    IOWrapper infile;
    infile.Open("test_microphysics.inp", IOWrapper::FileMode::read);

    pinput = new ParameterInput;
    pinput->LoadFromFile(infile);
    infile.Close();

    IndexMap::InitFromAthenaInput(pinput);
    Thermodynamics::InitFromAthenaInput(pinput);

    // set up mesh
    int restart = false;
    int mesh_only = false;
    pmesh = new Mesh(pinput, mesh_only);

    // set up components
    for (int b = 0; b < pmesh->nblocal; ++b) {
      MeshBlock *pmb = pmesh->my_blocks(b);
      pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pinput);
    }

    pmesh->Initialize(restart, pinput);
  }

  virtual void TearDown() {
    delete pinput;
    delete pmesh;
    Thermodynamics::Destroy();
    IndexMap::Destroy();
  }
};

TEST_F(TestMicrophysics, microphysics) {
  auto pmicro = pmesh->my_blocks(0)->pimpl->pmicro;

  EXPECT_EQ(pmicro->GetNumSystems(), 2);
  EXPECT_EQ(pmicro->GetSystemName(0), "water-system");
  EXPECT_EQ(pmicro->GetSystemName(1), "ammonia-system");
}

TEST_F(TestMicrophysics, assemble_matrix) {
  auto pindex = IndexMap::GetInstance();
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<AirParcel> air_column(1);
  auto &air = air_column[0];
  air.SetType(AirParcel::Type::MoleFrac);
  air.SetZero();

  int iH2O = pindex->GetSpeciesId("vapor.H2O");
  int iH2Oc = pindex->GetSpeciesId("cloud.H2O(c)");
  int iH2Op = pindex->GetSpeciesId("cloud.H2O(p)");

  int iNH3 = pindex->GetSpeciesId("vapor.NH3");
  int iNH3c = pindex->GetSpeciesId("cloud.NH3(c)");
  int iNH3p = pindex->GetSpeciesId("cloud.NH3(p)");

  EXPECT_EQ(iH2Oc, 7);
  EXPECT_EQ(iNH3c, 8);
  EXPECT_EQ(iH2Op, 9);
  EXPECT_EQ(iNH3p, 10);

  air.w[IDN] = 160.;
  air.w[IPR] = 7e4;
  air.w[iH2O] = 0.02;
  air.w[iNH3] = 0.30;

  pthermo->SaturationAdjustment(air_column);

  air.w[iNH3p] = air.w[iNH3c];
  air.w[iNH3c] = 0.;
  air.w[iNH3] = 0.;

  auto pmicro = pmesh->my_blocks(0)->pimpl->pmicro;

  auto psys_h2o = pmicro->GetSystem(0);
  psys_h2o->AssembleReactionMatrix(air, 0.);

  Real const *rate_h2o = static_cast<Kessler94 *>(psys_h2o.get())->GetRatePtr();
  Real const *jacb_h2o =
      static_cast<Kessler94 *>(psys_h2o.get())->GetJacobianPtr();

  EXPECT_NEAR(rate_h2o[0], 0., 1.e-8);
  EXPECT_NEAR(rate_h2o[1], -0.0199808, 1.e-8);
  EXPECT_NEAR(rate_h2o[2], 0.0199808, 1.e-8);

  EXPECT_NEAR(jacb_h2o[Kessler94::Size + 1], -1., 1.e-8);
  EXPECT_NEAR(jacb_h2o[2 * Kessler94::Size - 1], 1., 1.e-8);

  auto psys_nh3 = pmicro->GetSystem(1);
  psys_nh3->AssembleReactionMatrix(air, 0.);

  Real const *rate_nh3 = static_cast<Kessler94 *>(psys_nh3.get())->GetRatePtr();
  Real const *jacb_nh3 =
      static_cast<Kessler94 *>(psys_nh3.get())->GetJacobianPtr();

  EXPECT_NEAR(rate_nh3[0], 0.000834315, 1.e-8);
  EXPECT_NEAR(rate_nh3[1], 0., 1.e-8);
  EXPECT_NEAR(rate_nh3[2], -0.000834315, 1.e-8);

  EXPECT_NEAR(jacb_nh3[0], -0.00274292, 1.e-8);
  EXPECT_NEAR(jacb_nh3[Kessler94::Size - 1], 0.00274292, 1.e-8);

  EXPECT_NEAR(jacb_nh3[Kessler94::Size + 1], -1., 1.e-8);
  EXPECT_NEAR(jacb_nh3[2 * Kessler94::Size - 1], 1., 1.e-8);

  EXPECT_NEAR(jacb_nh3[2 * Kessler94::Size], 0.30417037565, 1.e-8);
  EXPECT_NEAR(jacb_nh3[3 * Kessler94::Size - 1], -0.30417037564, 1.e-8);
}

TEST_F(TestMicrophysics, evolve_one_step) {
  auto pindex = IndexMap::GetInstance();
  auto pthermo = Thermodynamics::GetInstance();

  std::vector<AirParcel> air_column(1);
  auto &air = air_column[0];
  air.SetType(AirParcel::Type::MoleFrac);
  air.SetZero();

  int iH2O = pindex->GetSpeciesId("vapor.H2O");
  int iH2Oc = pindex->GetSpeciesId("cloud.H2O(c)");
  int iH2Op = pindex->GetSpeciesId("cloud.H2O(p)");

  int iNH3 = pindex->GetSpeciesId("vapor.NH3");
  int iNH3c = pindex->GetSpeciesId("cloud.NH3(c)");
  int iNH3p = pindex->GetSpeciesId("cloud.NH3(p)");

  air.w[IDN] = 160.;
  air.w[IPR] = 7e4;
  air.w[iH2O] = 0.02;
  air.w[iNH3] = 0.30;

  pthermo->SaturationAdjustment(air_column);

  air.w[iNH3p] = air.w[iNH3c];
  air.w[iNH3c] = 0.;
  air.w[iNH3] = 0.;

  auto pmicro = pmesh->my_blocks(0)->pimpl->pmicro;

  auto psys_h2o = pmicro->GetSystem(0);
  psys_h2o->AssembleReactionMatrix(air, 0.);
  psys_h2o->EvolveOneStep(&air, 0., 1.e5);

  EXPECT_NEAR(air.w[iH2O], 1.98769e-05, 1e-8);
  EXPECT_NEAR(air.w[iH2Oc], 0., 1e-8);
  EXPECT_NEAR(air.w[iH2Op], 0.0199801231, 1e-8);

  auto psys_nh3 = pmicro->GetSystem(1);
  psys_nh3->AssembleReactionMatrix(air, 0.);
  psys_nh3->EvolveOneStep(&air, 0., 1.e5);
}

int main(int argc, char *argv[]) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);
  auto app = Application::GetInstance();

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
