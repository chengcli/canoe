#include "coordinates.hpp"

#include <torch/torch.h>

namespace canoe {

CoordinatesImpl::CoordinatesImpl() : options(CoordinatesOptions()) { reset(); }

CoordinatesImpl::CoordinatesImpl(const CoordinatesOptions& options_)
    : options(options_) {
  reset();
}

void CoordinatesImpl::reset() {
  auto nc1 = options.nx1() > 1 ? options.nx1() + 2 * options.nghost() : 1;
  auto nc2 = options.nx2() > 1 ? options.nx2() + 2 * options.nghost() : 1;
  auto nc3 = options.nx3() > 1 ? options.nx3() + 2 * options.nghost() : 1;

  register_buffer("area1", torch::empty({1, nc3, nc2, nc1 + 1}));
  register_buffer("area2", torch::empty({1, nc3, nc2 + 1, nc1}));
  register_buffer("area3", torch::empty({1, nc3 + 1, nc2, nc1}));
  register_buffer("vol", torch::empty({1, nc3, nc2, nc1}));
}

int64_t CoordinatesImpl::ncells3() const { return vol().size(1); }

int64_t CoordinatesImpl::ncells2() const { return vol().size(2); }

int64_t CoordinatesImpl::ncells1() const { return vol().size(3); }

torch::TensorList CoordinatesImpl::areas() const {
  auto buf = named_buffers(/*recurse=*/false);
  return {buf["area1"], buf["area2"], buf["area3"]};
}

torch::Tensor CoordinatesImpl::areas(int64_t dim) const {
  auto buf = named_buffers(/*recurse=*/false);

  switch (dim) {
    case 3:
      return buf["area1"];
    case 2:
      return buf["area2"];
    case 1:
      return buf["area3"];
    default:
      throw std::invalid_argument("dim must be 1, 2, or 3");
  }
}

torch::Tensor CoordinatesImpl::vol() const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["vol"];
}

torch::TensorList CoordinatesImpl::cos_theta() const {
  auto buf = named_buffers(/*recurse=*/false);
  return {buf["cos_theta1"], buf["cos_theta2"], buf["cos_theta3"]};
}

}  // namespace canoe
