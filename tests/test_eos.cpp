// external
#include <gtest/gtest.h>

// torch
#include <torch/torch.h>

// snap
#include <snap/eos/eos.hpp>

TEST(cons2prim_hydro_ideal, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  int64_t NHYDRO = 14;
  int64_t ncloud = 5;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto cons = torch::randn({NHYDRO, 1, 5, 5}, torch::device(torch::kMPS));

  auto gammad = torch::randn({1, 5, 5}, torch::device(torch::kMPS));
  auto rmu = torch::randn({nvapor + ncloud}, torch::device(torch::kMPS));
  auto rcv = torch::randn({nvapor + ncloud}, torch::device(torch::kMPS));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto prim = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);
  auto cons2 = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);

  EXPECT_TRUE(torch::allclose(cons, cons2, 1.E-6));
}

TEST(prim2cons_hydro_ideal, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  int64_t NHYDRO = 14;
  int64_t ncloud = 5;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto prim = torch::randn({NHYDRO, 1, 5, 5}, torch::device(torch::kMPS));

  auto gammad = torch::randn({1, 5, 5}, torch::device(torch::kMPS));
  auto rmu = torch::randn({nvapor + ncloud}, torch::device(torch::kMPS));
  auto rcv = torch::randn({nvapor + ncloud}, torch::device(torch::kMPS));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto cons = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);
  auto prim2 = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);

  EXPECT_TRUE(torch::allclose(prim, prim2, 1.E-6));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
