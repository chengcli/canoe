// C/C++
#include <algorithm>
#include <cmath>
#include <vector>

// external
#include <gtest/gtest.h>
#include <application/application.hpp>

// harp
#include <harp/read_stellar_flux.hpp>

std::string data_folder = "ck_data_01242024/";

TEST(read_stellar_flux, test_case1) {
  std::vector<double> expected_result = {
      1059434.8118773764, 6035941.569908597, 7520131.630087908,
      8488784.244992392,  5729875.721026806, 1400979.114668621,
      1083717.0198599927, 331262.2156410236, 296581.7003886011,
      44557.440115296864, 4298.3275481968085};

  auto app = Application::GetInstance();
  auto file1 = app->FindResource(data_folder + "sw_band_flux_HD189_11.txt");
  auto file2 = app->FindResource(data_folder + "wavelengths_GCM_11.txt");
  auto total_out = read_stellar_flux(file1, file2);

  EXPECT_EQ(total_out.first, expected_result);
}

TEST(read_stellar_flux, test_case2) {
  std::vector<double> expected_result = {0.260, 0.420, 0.610, 0.850,
                                         1.320, 2.020, 2.500, 3.500,
                                         4.400, 8.70,  20.00, 324.68};

  auto app = Application::GetInstance();
  auto file1 = app->FindResource(data_folder + "sw_band_flux_HD189_11.txt");
  auto file2 = app->FindResource(data_folder + "wavelengths_GCM_11.txt");
  auto total_out = read_stellar_flux(file1, file2);

  EXPECT_EQ(total_out.second, expected_result);
}

TEST(read_stellar_flux, test_case3) {
  std::vector<double> output1 = {
      1059434.8118773764, 6035941.569908597, 7520131.630087908,
      8488784.244992392,  5729875.721026806, 1400979.114668621,
      1083717.0198599927, 331262.2156410236, 296581.7003886011,
      44557.440115296864, 4298.3275481968085};

  std::vector<double> output2 = {0.260, 0.420, 0.610, 0.850, 1.320, 2.020,
                                 2.500, 3.500, 4.400, 8.70,  20.00, 324.68};

  auto expected_result = std::make_pair(output1, output2);

  auto app = Application::GetInstance();
  auto file1 = app->FindResource(data_folder + "sw_band_flux_HD189_11.txt");
  auto file2 = app->FindResource(data_folder + "wavelengths_GCM_11.txt");

  auto result = read_stellar_flux(file1, file2);
  EXPECT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
