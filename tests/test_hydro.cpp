// external
#include <gtest/gtest.h>

// torch
#include <torch/torch.h>

// snap
#include <snap/hydro/hydro.hpp>

using namespace canoe;

TEST(flux_divergence, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  std::vector<torch::Tensor> flux(3);
  std::vector<torch::Tensor> area(3);

  int64_t NHYDRO = 14;
  auto vol = torch::randn({1, 1, 5, 5},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));
  vol = vol.abs();

  for (int i = 0; i < 3; ++i) {
    flux[i] = torch::randn({NHYDRO, 1, 5, 5},
                           torch::device(torch::kMPS).dtype(torch::kFloat32));
    area[i] = torch::randn({1, 1, 5, 5},
                           torch::device(torch::kMPS).dtype(torch::kFloat32));
    area[i] = area[i].abs();
  }

  auto out = flux_divergence(flux, area, vol);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
