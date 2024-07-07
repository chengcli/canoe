// external
#include <gtest/gtest.h>

// torch
#include <torch/torch.h>

// snap
#include <snap/eos/eos.hpp>

using namespace canoe;

TEST(cons2prim_hydro_ideal, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  int64_t NHYDRO = 5;
  int64_t ncloud = 0;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto cons = torch::randn({NHYDRO, 1, 20, 20},
                           torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto gammad = torch::randn({1, 20, 20},
                             torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rmu = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rcv = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto start = std::chrono::high_resolution_clock::now();

  auto prim = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);
  auto cons2 = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time taken by test body: " << elapsed.count() << " seconds"
            << std::endl;

  std::cout << (cons - cons2).min() << std::endl;
  std::cout << (cons - cons2).max() << std::endl;

  EXPECT_TRUE(torch::allclose(cons, cons2, 1.E-4, 1.E-4));
}

TEST(prim2cons_hydro_ideal, mps_case1) {
  if (!torch::mps::is_available()) {
    GTEST_SKIP() << "MPS device is not available. Skipping test.";
  }

  int64_t NHYDRO = 14;
  int64_t ncloud = 5;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto prim = torch::randn({NHYDRO, 1, 5, 5},
                           torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto gammad = torch::randn({1, 5, 5},
                             torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rmu = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  auto rcv = torch::randn({nvapor + ncloud},
                          torch::device(torch::kMPS).dtype(torch::kFloat32));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto cons = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);
  auto prim2 = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);

  std::cout << (prim - prim2).abs().max() << std::endl;
  EXPECT_TRUE(torch::allclose(prim, prim2, 1.E-4, 1.E-4));
}

TEST(cons2prim_hydro_ideal, cpu_case1) {
  int64_t NHYDRO = 5;
  int64_t ncloud = 0;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto cons = torch::randn({NHYDRO, 1, 20, 20},
                           torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto gammad = torch::randn({1, 20, 20},
                             torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto rmu = torch::randn({nvapor + ncloud},
                          torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto rcv = torch::randn({nvapor + ncloud},
                          torch::device(torch::kCPU).dtype(torch::kFloat64));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto start = std::chrono::high_resolution_clock::now();

  auto prim = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);
  auto cons2 = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Time taken by test body: " << elapsed.count() << " seconds"
            << std::endl;

  std::cout << (cons - cons2).min() << std::endl;
  std::cout << (cons - cons2).max() << std::endl;

  EXPECT_TRUE(torch::allclose(cons, cons2, 1.E-5, 1.E-5));
}

TEST(cons2prim_hydro_ideal, cpu_case2) {
  int64_t NHYDRO = 5;
  int64_t ncloud = 0;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto cons = torch::randn({NHYDRO, 1, 20, 20},
                           torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto gammad = torch::randn({1, 20, 20},
                             torch::device(torch::kCPU).dtype(torch::kFloat64));

  gammad.normal_(0, 1);

  auto prim = eos_cons2prim_hydro_ideal(cons, gammad);
  auto cons2 = eos_prim2cons_hydro_ideal(prim, gammad);

  EXPECT_TRUE(torch::allclose(cons, cons2, 1.E-5, 1.E-5));
}

TEST(prim2cons_hydro_ideal, cpu_case1) {
  int64_t NHYDRO = 14;
  int64_t ncloud = 5;
  int64_t nvapor = NHYDRO - 5 - ncloud;

  auto prim = torch::randn({NHYDRO, 1, 5, 5},
                           torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto gammad = torch::randn({1, 5, 5},
                             torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto rmu = torch::randn({nvapor + ncloud},
                          torch::device(torch::kCPU).dtype(torch::kFloat64));

  auto rcv = torch::randn({nvapor + ncloud},
                          torch::device(torch::kCPU).dtype(torch::kFloat64));

  gammad.normal_(0, 1);
  rmu.normal_(0, 1);
  rcv.normal_(0, 1);

  auto cons = eos_prim2cons_hydro_ideal(prim, gammad, rmu, rcv, ncloud);
  auto prim2 = eos_cons2prim_hydro_ideal(cons, gammad, rmu, rcv, ncloud);

  std::cout << (prim - prim2).abs().max() << std::endl;
  EXPECT_TRUE(torch::allclose(prim, prim2, 1.E-5, 1.E-5));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
