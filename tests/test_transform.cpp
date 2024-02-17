// C/C++
#include <algorithm>
#include <cmath>

// external
#include <gtest/gtest.h>

// athena
#include <athena/athena.hpp>

// exo3
#include <exo3/cubed_sphere_utility.hpp>

namespace cs = CubedSphereUtility;

TEST(CovariantToContravariant, test1) {
  Real w[5] = {1., 2., 3., 4., 5.};
  Real cth = 0.5;

  cs::CovariantToContravariant(w, cth);

  EXPECT_DOUBLE_EQ(w[0], 1.);
  EXPECT_DOUBLE_EQ(w[1], 2.);
  EXPECT_DOUBLE_EQ(w[2], 4. / 3.);
  EXPECT_DOUBLE_EQ(w[3], 10. / 3.);
  EXPECT_DOUBLE_EQ(w[4], 5.);
}

TEST(ContravariantToOrthogonal, test1) {
  Real w[5] = {1., 2., 3., 4., 5.};
  Real y[5] = {1., 2., 3., 4., 5.};
  Real cth = 0.5;
  Real sth = sqrt(1. - cth * cth);

  cs::ContravariantToCovariant(y, cth);

  Real ek1 = w[IVX] * y[IVX] + w[IVY] * y[IVY] + w[IVZ] * y[IVZ];

  w[IVY] += w[IVZ] * cth;
  w[IVZ] *= sth;

  Real ek2 = w[IVX] * w[IVX] + w[IVY] * w[IVY] + w[IVZ] * w[IVZ];

  EXPECT_DOUBLE_EQ(w[0], 1.);
  EXPECT_DOUBLE_EQ(w[1], 2.);
  EXPECT_DOUBLE_EQ(w[2], 5);
  EXPECT_DOUBLE_EQ(w[3], sqrt(12));
  EXPECT_DOUBLE_EQ(w[4], 5.);
  EXPECT_DOUBLE_EQ(ek1, ek2);

  w[IVZ] /= sth;
  w[IVY] -= w[IVZ] * cth;

  EXPECT_DOUBLE_EQ(w[0], 1.);
  EXPECT_DOUBLE_EQ(w[1], 2.);
  EXPECT_DOUBLE_EQ(w[2], 3.);
  EXPECT_DOUBLE_EQ(w[3], 4.);
  EXPECT_DOUBLE_EQ(w[4], 5.);
}

TEST(ContravariantToCovariant, test1) {
  Real w[5] = {1., 2., 4. / 3., 10. / 3., 5.};
  Real cth = 0.5;

  cs::ContravariantToCovariant(w, cth);

  EXPECT_DOUBLE_EQ(w[0], 1.);
  EXPECT_DOUBLE_EQ(w[1], 2.);
  EXPECT_DOUBLE_EQ(w[2], 3.);
  EXPECT_DOUBLE_EQ(w[3], 4.);
  EXPECT_DOUBLE_EQ(w[4], 5.);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

  return result;
}
