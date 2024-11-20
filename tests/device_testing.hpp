#pragma once

// external
#include <gtest/gtest.h>

// torch
#include <torch/torch.h>

struct Parameters {
  torch::DeviceType device_type;
  torch::Dtype dtype;
};

void PrintTo(const Parameters& param, std::ostream* os) {
  std::string device_str = torch::Device(param.device_type).str();
  std::string dtype_str = torch::toString(param.dtype);
  *os << "Device: " << device_str << ", Dtype: " << dtype_str;
}

class DeviceTest : public testing::TestWithParam<Parameters> {
 protected:
  torch::Device device = torch::kCPU;
  torch::Dtype dtype = torch::kFloat32;

  void SetUp() override {
    // Get the current parameters
    auto param = GetParam();
    device = torch::Device(param.device_type);
    dtype = param.dtype;

    // Check if the device is available, and skip the test if not
    if (device.type() == torch::kCUDA && !torch::cuda::is_available()) {
      GTEST_SKIP() << "CUDA is not available, skipping test.";
    }

    if (device.type() == torch::kMPS && !torch::hasMPS()) {
      GTEST_SKIP() << "MPS is not available, skipping test.";
    }
  }
};

INSTANTIATE_TEST_SUITE_P(
    DeviceAndDtype, DeviceTest,
    testing::Values(Parameters{torch::kCPU, torch::kFloat32},
                    Parameters{torch::kCPU, torch::kFloat64},
                    Parameters{torch::kMPS, torch::kFloat32},
                    Parameters{torch::kCUDA, torch::kFloat32},
                    Parameters{torch::kCUDA, torch::kFloat64}),
    [](const testing::TestParamInfo<DeviceTest::ParamType>& info) {
      std::string name = torch::Device(info.param.device_type).str();
      name += "_";
      name += torch::toString(info.param.dtype);
      std::replace(name.begin(), name.end(), '.', '_');
      return name;
    });
