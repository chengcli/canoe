// C/C++
#include <algorithm>
#include <cmath>

// external
#include <gtest/gtest.h>

// athena
#include <athena/reconstruct/interpolation.hpp>

// torch
#include <torch/torch.h>

// snap
#include <snap/reconstruct/interpolation.hpp>

using namespace canoe;

TEST(interp_weno3, test_case1) {
  double phim1 = 1.0;
  double phi = 2.0;
  double phip1 = 3.0;
  double result = interp_weno3(phim1, phi, phip1);
  double expected_result = 1.5;
  EXPECT_NEAR(result, expected_result, 1.E-10);
}

TEST(interp_weno5, test_case1) {
  double phim2 = 1.0;
  double phim1 = 2.0;
  double phi = 3.0;
  double phip1 = 4.0;
  double phip2 = 5.0;
  double result = interp_weno5(phim2, phim1, phi, phip1, phip2);
  double expected_result = 2.5000000000000004;
  EXPECT_NEAR(result, expected_result, 1.E-10);
}

TEST(interp_weno5m_torch, test_case2) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({5});
    float result1 = interp_weno5(phi[0].item<float>(), phi[1].item<float>(),
                                 phi[2].item<float>(), phi[3].item<float>(),
                                 phi[4].item<float>());
    float result2 = interp.left(phi).item<float>();
    EXPECT_NEAR(result1, result2, 2.E-6);
  }
}

TEST(interp_weno5m_torch, test_case3) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({2, 5});
    float result1 =
        interp_weno5(phi[0][0].item<float>(), phi[0][1].item<float>(),
                     phi[0][2].item<float>(), phi[0][3].item<float>(),
                     phi[0][4].item<float>());
    float result2 =
        interp_weno5(phi[1][0].item<float>(), phi[1][1].item<float>(),
                     phi[1][2].item<float>(), phi[1][3].item<float>(),
                     phi[1][4].item<float>());
    torch::Tensor result = interp.left(phi);
    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

TEST(interp_weno5m_torch, test_case4) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({5, 2});
    float result1 =
        interp_weno5(phi[0][0].item<float>(), phi[1][0].item<float>(),
                     phi[2][0].item<float>(), phi[3][0].item<float>(),
                     phi[4][0].item<float>());
    float result2 =
        interp_weno5(phi[0][1].item<float>(), phi[1][1].item<float>(),
                     phi[2][1].item<float>(), phi[3][1].item<float>(),
                     phi[4][1].item<float>());
    auto phiu = phi.unfold(0, 5, 1);
    torch::Tensor result = interp.left(phiu).squeeze(0);

    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

TEST(interp_weno5p_torch, test_case5) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({5});
    float result1 = interp_weno5(phi[4].item<float>(), phi[3].item<float>(),
                                 phi[2].item<float>(), phi[1].item<float>(),
                                 phi[0].item<float>());
    float result2 = interp.right(phi).item<float>();
    EXPECT_NEAR(result1, result2, 2.E-6);
  }
}

TEST(interp_weno5p_torch, test_case6) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({2, 5});
    float result1 =
        interp_weno5(phi[0][4].item<float>(), phi[0][3].item<float>(),
                     phi[0][2].item<float>(), phi[0][1].item<float>(),
                     phi[0][0].item<float>());
    float result2 =
        interp_weno5(phi[1][4].item<float>(), phi[1][3].item<float>(),
                     phi[1][2].item<float>(), phi[1][1].item<float>(),
                     phi[1][0].item<float>());
    torch::Tensor result = interp.right(phi);
    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

TEST(interp_weno5p_torch, test_case7) {
  Weno5Interp interp;
  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({5, 2});
    float result1 =
        interp_weno5(phi[4][0].item<float>(), phi[3][0].item<float>(),
                     phi[2][0].item<float>(), phi[1][0].item<float>(),
                     phi[0][0].item<float>());
    float result2 =
        interp_weno5(phi[4][1].item<float>(), phi[3][1].item<float>(),
                     phi[2][1].item<float>(), phi[1][1].item<float>(),
                     phi[0][1].item<float>());
    auto phiu = phi.unfold(0, 5, 1);
    torch::Tensor result = interp.right(phiu).squeeze(0);
    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

TEST(interp_weno5m_torch, test_case_mps1) {
  if (!torch::mps::is_available()) {
    return;
  }
  Weno5Interp interp(torch::kMPS);

  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({2, 5}, torch::device(torch::kMPS));

    float result1 =
        interp_weno5(phi[0][0].item<float>(), phi[0][1].item<float>(),
                     phi[0][2].item<float>(), phi[0][3].item<float>(),
                     phi[0][4].item<float>());

    float result2 =
        interp_weno5(phi[1][0].item<float>(), phi[1][1].item<float>(),
                     phi[1][2].item<float>(), phi[1][3].item<float>(),
                     phi[1][4].item<float>());

    torch::Tensor result = interp.left(phi);
    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

TEST(interp_weno5p_torch, test_case_mps2) {
  Weno5Interp interp(torch::kMPS);

  for (int i = 0; i < 10; ++i) {
    torch::Tensor phi = torch::randn({5, 2}, torch::device(torch::kMPS));

    float result1 =
        interp_weno5(phi[4][0].item<float>(), phi[3][0].item<float>(),
                     phi[2][0].item<float>(), phi[1][0].item<float>(),
                     phi[0][0].item<float>());
    float result2 =
        interp_weno5(phi[4][1].item<float>(), phi[3][1].item<float>(),
                     phi[2][1].item<float>(), phi[1][1].item<float>(),
                     phi[0][1].item<float>());

    auto phiu = phi.unfold(0, 5, 1);
    torch::Tensor result = interp.right(phiu).squeeze(0);
    EXPECT_NEAR(result1, result[0].item<float>(), 2.E-6);
    EXPECT_NEAR(result2, result[1].item<float>(), 2.E-6);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
