// C/C++
#include <algorithm>
#include <cmath>
#include <vector>

// external
#include <gtest/gtest.h>

#include <application/application.hpp>

// athena
#include <athena/athena_arrays.hpp>

// exo3
#include <exo3/cs_velocity_rotation.hpp>

TEST(vel_zab_from_p1_test, test_case_to_4) {
  std::vector <Real> result = {1,2,3};
  std::vector <Real> expected_result = {-2,-3,1};
  CubedSphereUtility::vel_zab_from_p1_test(result, 4);
  EXPECT_EQ(result, expected_result);
}

TEST(vel_zab_from_p1_test, test_case_to_6) {
  std::vector <Real> result = {1,2,3};
  std::vector <Real> expected_result = {-3,2,1};
  CubedSphereUtility::vel_zab_from_p1_test(result, 6);
  EXPECT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}