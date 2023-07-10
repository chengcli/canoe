// external
#include <gtest/gtest.h>

// application
#include <application/application.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/field/field.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/scalars/scalars.hpp>

// canoe
#include <impl.hpp>
#include <index_map.hpp>
#include <variable.hpp>

class TestImpl : public testing::Test {
 protected:
  ParameterInput *pinput;
  Mesh *pmesh;

  virtual void SetUp() {
    // code here will execute just before the test ensues
    IOWrapper infile;
    infile.Open("test_thermodynamics.inp", IOWrapper::FileMode::read);

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

    // set up variables
    auto pmb = pmesh->my_blocks(0);
    auto phydro = pmb->phydro;
    auto pcloud = pmb->pimpl->pcloud;
    auto ptracer = pmb->pimpl->ptracer;
    auto pfield = pmb->pfield;
    auto pcoord = pmb->pcoord;
    auto pscalars = pmb->pscalars;
    auto peos = pmb->peos;

    int ks = pmb->ks, js = pmb->js, is = pmb->is;

    phydro->w(IDN, ks, js, is) = 1.;
    phydro->w(IPR, ks, js, is) = 10.;
    phydro->w(IVX, ks, js, is) = 0.1;
    phydro->w(IVY, ks, js, is) = 0.2;
    phydro->w(IVZ, ks, js, is) = 0.3;

    for (int n = 1; n <= NVAPOR; ++n)
      phydro->w(n, ks, js, is) =
          1. / (1 + NVAPOR + 2 * NCLOUD + 3 * NCHEMISTRY);

    for (int n = 0; n < NCLOUD; ++n)
      pcloud->w(n, ks, js, is) =
          2. / (1 + NVAPOR + 2 * NCLOUD + 3 * NCHEMISTRY);

    for (int n = 0; n < NTRACER; ++n) ptracer->w(n, ks, js, is) = 2.;

    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is,
                               is, js, js, ks, ks);

    peos->PassiveScalarPrimitiveToConserved(pscalars->r, phydro->u, pscalars->s,
                                            pcoord, is, is, js, js, ks, ks);
  }

  virtual void TearDown() {
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be

    delete pinput;
    delete pmesh;
    Thermodynamics::Destroy();
    IndexMap::Destroy();
  }
};

TEST_F(TestImpl, GatherPrimitive) {
  Variable var(Variable::Type::MoleFrac);
  auto pmb = pmesh->my_blocks(0);
  auto pimpl = pmb->pimpl;
  auto phydro = pmb->phydro;
  auto pcloud = pimpl->pcloud;
  auto ptracer = pimpl->ptracer;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  pimpl->GatherFromPrimitive(&var, ks, js, is);
  pimpl->DistributeToPrimitive(var, ks, js, is);

  EXPECT_DOUBLE_EQ(phydro->w(IDN, ks, js, is), 1.);
  EXPECT_DOUBLE_EQ(phydro->w(IPR, ks, js, is), 10.);
  EXPECT_DOUBLE_EQ(phydro->w(IVX, ks, js, is), 0.1);
  EXPECT_DOUBLE_EQ(phydro->w(IVY, ks, js, is), 0.2);
  EXPECT_DOUBLE_EQ(phydro->w(IVZ, ks, js, is), 0.3);

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_DOUBLE_EQ(phydro->w(n, ks, js, is),
                     1. / (1 + NVAPOR + 2 * NCLOUD + 3 * NCHEMISTRY));

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_DOUBLE_EQ(pcloud->w(n, ks, js, is),
                     2. / (1 + NVAPOR + 2 * NCLOUD + 3 * NCHEMISTRY));

  for (int n = 0; n < NTRACER; ++n)
    EXPECT_DOUBLE_EQ(ptracer->w(n, ks, js, is), 2.);
}

TEST_F(TestImpl, GatherConserved) {
  Variable var(Variable::Type::MassConc);
  auto pmb = pmesh->my_blocks(0);
  auto pimpl = pmb->pimpl;
  auto phydro = pmb->phydro;
  auto pcloud = pimpl->pcloud;
  auto ptracer = pimpl->ptracer;

  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  pimpl->GatherFromConserved(&var, ks, js, is);
  pimpl->DistributeToConserved(var, ks, js, is);

  EXPECT_NEAR(phydro->u(IDN, ks, js, is), 0.818181818181, 1e-10);
  EXPECT_NEAR(phydro->u(IEN, ks, js, is), 24.98002735832, 1e-8);

  EXPECT_DOUBLE_EQ(phydro->u(IVX, ks, js, is), 0.1);
  EXPECT_DOUBLE_EQ(phydro->u(IVY, ks, js, is), 0.2);
  EXPECT_DOUBLE_EQ(phydro->u(IVZ, ks, js, is), 0.3);

  for (int n = 1; n <= NVAPOR; ++n)
    EXPECT_NEAR(phydro->u(n, ks, js, is), 0.0909091, 1e-8);

  for (int n = 0; n < NCLOUD; ++n)
    EXPECT_NEAR(pcloud->u(n, ks, js, is), 0.1487603305785, 1e-10);

  for (int n = 0; n < NTRACER; ++n)
    EXPECT_NEAR(ptracer->u(n, ks, js, is), 1.6363636363636, 1e-10);
};

int main(int argc, char **argv) {
  Application::Start(argc, argv);

  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  Application::Destroy();

  return result;
}
