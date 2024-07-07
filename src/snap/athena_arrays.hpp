#pragma once
// C/C++
#include <exception>

// athena
#include <athena/athena_arrays.hpp>

// torch
#include <torch/torch.h>

template <typename T>
void AthenaArray<T>::toDevice(c10::DeviceType device) {
  if (nx5_ != 1 || nx6_ != 1) {
    throw std::runtime_error("AthenaArray::toDevice: nx5 and nx6 must be 1");
  }

  int64_t str1 = 1;
  int64_t str2 = nx1_;
  int64_t str3 = nx2_ * nx1_;
  int64_t str4 = nx3_ * nx2_ * nx1_;

  ptensor_ = std::make_shared<torch::Tensor>(
      torch::from_blob(pdata_, {nx4_, nx3_, nx2_, nx1_},
                       {str4, str3, str2, str1}, nullptr,
                       torch::dtype(torch::kFloat32))
          .to(device));
}

template <typename T>
void AthenaArray<T>::fromDevice() {
  int64_t str1 = 1;
  int64_t str2 = nx1_;
  int64_t str3 = nx2_ * nx1_;
  int64_t str4 = nx3_ * nx2_ * nx1_;

  // create a temporary tensor holder
  torch::Tensor tmp = torch::from_blob(pdata_, {nx4_, nx3_, nx2_, nx1_},
                                       {str4, str3, str2, str1}, nullptr,
                                       torch::dtype(torch::kFloat32));

  tmp.copy_(*ptensor_);
}

template <typename T>
torch::Tensor& AthenaArray<T>::tensor() {
  return *ptensor_;
}

template <typename T>
torch::Tensor const& AthenaArray<T>::tensor() const {
  return *ptensor_;
}
