// external
#include <gtest/gtest.h>

// torch
#include <torch/torch.h>

// snap
#include <snap/riemann/riemann.hpp>

enum {
  DIM1 = 3,
  DIM2 = 2,
  DIM3 = 1,
};

using namespace canoe;

TEST(hydro_lmars, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  int64_t NHYDRO = 5;
  int64_t ncloud = 0;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto wl = torch::randn({NHYDRO, 1, 20, 20},
                         torch::device(torch::kMPS).dtype(torch::kFloat32));
  wl = wl.abs();
  auto wr = torch::randn({NHYDRO, 1, 20, 20},
                         torch::device(torch::kMPS).dtype(torch::kFloat32));
  wr = wr.abs();

  auto gammad = torch::randn({1, 20, 20},
                             torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rmu = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rcv = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  gammad = gammad.normal_(0, 1).abs();
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto start = std::chrono::high_resolution_clock::now();

  auto flux = rs_hydro_lmars(DIM1, wl, wr, gammad, rmu, rcv, ncloud);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time taken by test body: " << elapsed.count() << " seconds"
            << std::endl;
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
