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

TEST(vel_zab_to_zxy, test_ab_to_xy_to_ab) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zxy_to_zab(vz, vx, vy, PI/3, PI/4);
  std::cout << *vz << " " << *vx << " " << *vy;
  CubedSphereUtility::vel_zab_to_zxy(vz, vx, vy, PI/3, PI/4);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_1_to_2_to_1) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 2);
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 1);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_1_to_3_to_1) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 3);
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 1);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_1_to_4_to_1) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 4);
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 1);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_1_to_6_to_1) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 6);
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 1);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_2_to_1_to_2) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 1);
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 2);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_2_to_3_to_2) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 3);
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 2);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_2_to_4_to_2) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 4);
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 2);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_2_to_5_to_2) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 5);
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 2);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_3_to_1_to_3) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 1);
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_3_to_2_to_3) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 2);
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_3_to_5_to_3) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 5);
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_3_to_6_to_3) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 6);
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_4_to_2_to_4) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 2);
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 4);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_4_to_1_to_4) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 1);
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 4);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_4_to_5_to_4) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 5);
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 4);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_4_to_6_to_4) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 6);
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 4);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_5_to_2_to_5) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 2);
  CubedSphereUtility::vel_zab_from_p2(vz, vx, vy, PI/5*2, PI/8*3, 5);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_5_to_3_to_5) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 3);
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 5);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_5_to_4_to_5) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 4);
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 5);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_5_to_6_to_5) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 6);
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 5);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_6_to_1_to_6) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 1);
  CubedSphereUtility::vel_zab_from_p1(vz, vx, vy, PI/5*2, PI/8*3, 6);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_6_to_3_to_6) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 3);
  CubedSphereUtility::vel_zab_from_p3(vz, vx, vy, PI/5*2, PI/8*3, 6);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_6_to_4_to_6) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 4);
  CubedSphereUtility::vel_zab_from_p4(vz, vx, vy, PI/5*2, PI/8*3, 6);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

TEST(vel_zab_from_test, test_case_6_to_5_to_6) {
  Real result[3] = {1,2,3};
  Real *vz = result;
  Real *vx = result + 1;
  Real *vy = result + 2;
  Real expected_result[3] = {1,2,3};
  CubedSphereUtility::vel_zab_from_p6(vz, vx, vy, PI/5*2, PI/8*3, 5);
  CubedSphereUtility::vel_zab_from_p5(vz, vx, vy, PI/5*2, PI/8*3, 6);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(result[i], expected_result[i], 1e-5);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}