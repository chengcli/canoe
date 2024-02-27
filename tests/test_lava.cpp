// C/C++
#include <algorithm>
#include <cmath>

// external
#include <gtest/gtest.h>

// athena
#include <athena/reconstruct/interpolation.hpp>

// test temp, theta, thetav, mse of each air parcel

//expected results
double *temp_ptr = atm["TEMP"].data();
double *theta_ptr = atm["THETA"].data();
double *thetav_ptr = atm["THETAV"].data();
double *mse_ptr = atm["MSE"].data();

TEST(interp_weno3, test_case1) {
  int atm_size = 80;
  for (int i = 0; i < atm_size; ++i) {
    // results
    double temp = ;
    double theta = ;
    double thetav = ;
    double mse = ;
    
    EXPECT_NEAR(temp, temp_ptr[i], 1.E-10);
    EXPECT_NEAR(theta, theta_ptr[i], 1.E-10);
    EXPECT_NEAR(thetav, thetav_ptr[i], 1.E-10);
    EXPECT_NEAR(mse, mse_ptr[i], 1.E-10);
  }
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
