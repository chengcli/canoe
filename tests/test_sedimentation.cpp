// external
#include <gtest/gtest.h>

// athena
#include <athena/athena.hpp>

// microphysics
#include <microphysics/sedimentation.hpp>

// tests
#include "device_testing.hpp"

TEST_P(DeviceTest, forward) {
  auto options = SedimentationOptions();
  auto sed = Sedimentation(options);

  int nhydro = 5;
  int nx3 = 1;
  int nx2 = 1;
  int nx1 = 9;

  sed->to(device, dtype);

  auto w =
      torch::randn({nhydro, nx3, nx2, nx1}, torch::device(device).dtype(dtype));
  w[IDN] = w[IDN].abs();
  w[IPR] = w[IPR].abs();

  auto diag = std::make_shared<SharedData::element_type>();
  (*diag)["temperature"] =
      std::async(std::launch::async, [&]() {
        return 300. *
               torch::randn({nx3, nx2, nx1}, torch::device(device).dtype(dtype))
                   .abs();
      }).share();
  sed->set_shared_data(diag);

  auto vel = sed->forward(w);
  std::cout << "vel = " << vel << std::endl;
}

TEST_P(DeviceTest, two_species_forward) {
  auto options = SedimentationOptions();
  options.radius({0.0e-6, 10.e-6});
  options.density({1.0e3, 1.0e3});

  auto sed = Sedimentation(options);

  int nhydro = 5;
  int nx3 = 1;
  int nx2 = 1;
  int nx1 = 9;

  sed->to(device, dtype);

  auto w =
      torch::randn({nhydro, nx3, nx2, nx1}, torch::device(device).dtype(dtype));
  w[IDN] = w[IDN].abs();
  w[IPR] = w[IPR].abs();

  auto diag = std::make_shared<SharedData::element_type>();
  (*diag)["temperature"] =
      std::async(std::launch::async, [&]() {
        return 300. *
               torch::randn({nx3, nx2, nx1}, torch::device(device).dtype(dtype))
                   .abs();
      }).share();
  sed->set_shared_data(diag);

  auto vel = sed->forward(w);
  std::cout << "vel = " << vel << std::endl;
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // start_logging(argc, argv);

  return RUN_ALL_TESTS();
}
