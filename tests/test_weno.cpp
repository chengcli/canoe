// athena
#include <athena/reconstruct/interpolation.hpp>

// external
#include <gtest/gtest.h>

#include <cmath>

TEST(interp_weno3, test_case1) {
  double phim1 = 1.0;
  double phi = 2.0;
  double phip1 = 3.0;
  double result = interp_weno3(phim1, phi, phip1);
  double expected_result = 1.5;
  EXPECT_NEAR(result, expected_result, 1.E-10);
}

TEST(interp_weno5, test_case1) {
  double phim2 = 1.0;
  double phim1 = 2.0;
  double phi = 3.0;
  double phip1 = 4.0;
  double phip2 = 5.0;
  double result = interp_weno5(phim2, phim1, phi, phip1, phip2);
  double expected_result = 2.5000000000000004;
  EXPECT_NEAR(result, expected_result, 1.E-10);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
