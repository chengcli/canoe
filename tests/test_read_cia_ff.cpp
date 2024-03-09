// C/C++
#include <algorithm>
#include <cmath>
#include <vector>

// external
#include <gtest/gtest.h>

#include <application/application.hpp>

// athena
#include <athena/athena_arrays.hpp>

// opacity
#include <opacity/read_cia_ff.hpp>

std::string data_folder = "ck_data_01242024/cia/";

TEST(read_cia_reform, test_case1) {
  std::vector<double> expected_result = {4.47996926605873e-45,
                                         1.868712588877408e-44,
                                         1.2799236924204295e-44,
                                         4.467454212770098e-45,
                                         1.2948567646502106e-45,
                                         2.802009368955379e-46,
                                         2.308620846855808e-47,
                                         3.3960726642220515e-49,
                                         1e-99,
                                         1e-99,
                                         1e-99};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "He-H_reform.cia");
  auto data_table = read_cia_reform(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 11; ++i) {
    result.push_back(data(2, i));
  }

  EXPECT_EQ(result, expected_result);
}

TEST(read_cia_reform, test_case2) {
  std::vector<double> expected_result = {2.0649974391687292e-45,
                                         1.504861456948198e-45,
                                         9.0902488843771e-47,
                                         5.624908034565523e-48,
                                         4.523257975006434e-47,
                                         3.004769764226695e-45,
                                         6.423470305957112e-48,
                                         2.647376247478312e-49,
                                         1.9496282953781007e-51,
                                         1e-99,
                                         1e-99};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-He_reform.cia");
  auto data_table = read_cia_reform(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 11; ++i) {
    result.push_back(data(8, i));
  }

  EXPECT_EQ(result, expected_result);
}

TEST(read_cia_reform, test_case3) {
  std::vector<double> expected_result = {
      200.0,  225.0,  250.0,  275.0,  300.0,  325.0,  350.0,  375.0,  400.0,
      425.0,  450.0,  475.0,  500.0,  525.0,  550.0,  575.0,  600.0,  625.0,
      650.0,  675.0,  700.0,  725.0,  750.0,  775.0,  800.0,  825.0,  850.0,
      875.0,  900.0,  925.0,  950.0,  975.0,  1000.0, 1025.0, 1050.0, 1075.0,
      1100.0, 1125.0, 1150.0, 1175.0, 1200.0, 1225.0, 1250.0, 1275.0, 1300.0,
      1325.0, 1350.0, 1375.0, 1400.0, 1425.0, 1450.0, 1475.0, 1500.0, 1525.0,
      1550.0, 1575.0, 1600.0, 1625.0, 1650.0, 1675.0, 1700.0, 1725.0, 1750.0,
      1775.0, 1800.0, 1825.0, 1850.0, 1875.0, 1900.0, 1925.0, 1950.0, 1975.0,
      2000.0, 2025.0, 2050.0, 2075.0, 2100.0, 2125.0, 2150.0, 2175.0, 2200.0,
      2225.0, 2250.0, 2275.0, 2300.0, 2325.0, 2350.0, 2375.0, 2400.0, 2425.0,
      2450.0, 2475.0, 2500.0, 2525.0, 2550.0, 2575.0, 2600.0, 2625.0, 2650.0,
      2675.0, 2700.0, 2725.0, 2750.0, 2775.0, 2800.0, 2825.0, 2850.0, 2875.0,
      2900.0, 2925.0, 2950.0, 2975.0, 3000.0, 3025.0, 3050.0, 3075.0, 3100.0,
      3125.0, 3150.0, 3175.0, 3200.0, 3225.0, 3250.0, 3275.0, 3300.0, 3325.0,
      3350.0, 3375.0, 3400.0, 3425.0, 3450.0, 3475.0, 3500.0, 3525.0, 3550.0,
      3575.0, 3600.0, 3625.0, 3650.0, 3675.0, 3700.0, 3725.0, 3750.0, 3775.0,
      3800.0, 3825.0, 3850.0, 3875.0, 3900.0, 3925.0, 3950.0, 3975.0, 4000.0,
      4025.0, 4050.0, 4075.0, 4100.0, 4125.0, 4150.0, 4175.0, 4200.0, 4225.0,
      4250.0, 4275.0, 4300.0, 4325.0, 4350.0, 4375.0, 4400.0, 4425.0, 4450.0,
      4475.0, 4500.0, 4525.0, 4550.0, 4575.0, 4600.0, 4625.0, 4650.0, 4675.0,
      4700.0, 4725.0, 4750.0, 4775.0, 4800.0, 4825.0, 4850.0, 4875.0, 4900.0,
      4925.0, 4950.0, 4975.0, 5000.0, 5025.0, 5050.0, 5075.0, 5100.0, 5125.0,
      5150.0, 5175.0, 5200.0, 5225.0, 5250.0, 5275.0, 5300.0, 5325.0, 5350.0,
      5375.0, 5400.0, 5425.0, 5450.0, 5475.0, 5500.0, 5550.0, 5575.0, 5600.0,
      5625.0, 5750.0, 5775.0, 5800.0, 5825.0, 5850.0, 5875.0, 5900.0, 5925.0,
      5950.0, 5975.0, 6000.0, 6025.0, 6050.0, 6075.0, 6100.0, 6125.0, 6150.0,
      6175.0, 6200.0, 6225.0, 6250.0, 6275.0, 6300.0, 6325.0, 6350.0, 6375.0,
      6400.0, 6425.0, 6450.0, 6475.0, 6500.0, 6525.0, 6550.0, 6575.0, 6600.0,
      6625.0, 6650.0, 6675.0, 6700.0, 6725.0, 6750.0, 6775.0, 6800.0, 6825.0,
      6850.0, 6875.0, 6900.0, 6925.0, 6950.0, 6975.0, 7000.0, 7025.0, 7050.0,
      7075.0, 7100.0, 7125.0, 7150.0, 7175.0, 7275.0, 7300.0, 7325.0, 7350.0,
      7375.0, 7400.0, 7425.0, 7450.0, 7475.0, 7500.0, 7525.0, 7550.0, 7575.0,
      7600.0, 7625.0, 7650.0, 7675.0, 7700.0, 7725.0, 7750.0, 7775.0, 7800.0,
      7825.0, 7850.0, 7875.0, 7900.0, 7925.0, 7950.0, 7975.0, 8000.0, 8050.0,
      8100.0, 8150.0, 8200.0, 8250.0, 8300.0, 8350.0, 8400.0, 8450.0, 8500.0,
      8550.0, 8600.0, 8650.0, 8700.0, 8750.0, 8800.0, 8850.0, 8900.0, 8950.0,
      9000.0, 9100.0, 9200.0, 9300.0, 9400.0, 9500.0, 9600.0, 9700.0, 9800.0,
      9900.0};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-He_reform.cia");
  auto data_table = read_cia_reform(file);
  auto temperature_data = std::get<1>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 334; ++i) {
    result.push_back(temperature_data[i]);
  }
  EXPECT_EQ(result, expected_result);
}

TEST(read_cia_reform, test_case4) {
  std::vector<double> expected_result = {
      265.3997782431933,  824.712643678161,   1711.0762800417974,
      2564.935064935065,  3428.5714285714284, 4475.247524752475,
      6263.126312631262,  9670.231729055258,  14079.07425265188,
      20101.483216237313, 31135.531135531135};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-He_reform.cia");
  auto data_table = read_cia_reform(file);
  auto spectral_data = std::get<2>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 11; ++i) {
    result.push_back(spectral_data[i]);
  }
  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case1) {
  std::vector<double> expected_result = {7.74e-2, 1.07e-1, 1.21e-1, 1.34e-1,
                                         1.57e-1, 1.81e-1, 2.29e-1, 2.79e-1};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "He-_ff.txt");
  auto data_table = read_freefree(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(4, i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case2) {
  std::vector<double> expected_result = {8.70e-2, 1.24e-1, 1.46e-1, 1.67e-1,
                                         2.10e-1, 2.53e-1, 3.39e-1, 4.27e-1};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-_ff.txt");
  auto data_table = read_freefree(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(2, i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case3) {
  std::vector<double> expected_result = {2.80e1, 3.62e1, 4.08e1, 4.49e1,
                                         5.26e1, 5.98e1, 7.27e1, 8.40e1};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "He-_ff.txt");
  auto data_table = read_freefree(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(15, i));
  }
  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case4) {
  std::vector<double> expected_result = {7.16e1, 9.23e1, 1.01e2, 1.08e2,
                                         1.18e2, 1.26e2, 1.38e2, 1.47e2};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-_ff.txt");
  auto data_table = read_freefree(file);
  auto data = std::get<0>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(data(17, i));
  }

  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case5) {
  std::vector<double> expected_result = {0.5, 0.8, 1.0, 1.2,
                                         1.6, 2.0, 2.8, 3.6};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-_ff.txt");
  auto data_table = read_freefree(file);
  auto temperature_data = std::get<1>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 8; ++i) {
    result.push_back(temperature_data[i]);
  }

  EXPECT_EQ(result, expected_result);
}

TEST(read_freefree, test_case6) {
  std::vector<double> expected_result = {
      3505,  4142,  5063,  5696,  6509,  7594,  9113,  11391,  15188,
      18226, 22783, 30377, 36452, 45565, 60753, 91130, 113913, 151883};

  auto app = Application::GetInstance();
  auto file = app->FindResource(data_folder + "H2-_ff.txt");
  auto data_table = read_freefree(file);
  auto spectral_data = std::get<2>(data_table);

  std::vector<double> result;
  for (int i = 0; i < 18; ++i) {
    result.push_back(spectral_data[i]);
  }

  EXPECT_EQ(result, expected_result);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
