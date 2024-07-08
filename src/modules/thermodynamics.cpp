#include "thermodynamics.hpp"

#include <torch/torch.h>

namespace canoe {

ThermodynamicsImpl::ThermodynamicsImpl() : options(ThermodynamicsOptions()) {
  reset();
}

ThermodynamicsImpl::ThermodynamicsImpl(const ThermodynamicsOptions& options_)
    : options(options_) {
  reset();
}

void ThermodynamicsImpl::reset() {
  // Do nothing
  register_buffer("cv_ratio_mass", torch::empty({1 + nvapor() + ncloud()}));
  register_buffer("mu_ratio", torch::empty({1 + nvapor() + ncloud()}));
}

int64_t ThermodynamicsImpl::nvapor() const { return options.nvapor(); }

int64_t ThermodynamicsImpl::ncloud() const { return options.ncloud(); }

torch::Tensor ThermodynamicsImpl::mu_ratio() const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["mu_ratio"];
}

torch::Tensor ThermodynamicsImpl::cv_ratio_mass() const {
  auto buf = named_buffers(/*recurse=*/false);
  return buf["cv_ratio_mass"];
}

torch::Tensor ThermodynamicsImpl::gammad(torch::Tensor const& hydro_u) const {
  return torch::ones_like(hydro_u) * options.gammad_ref();
}

}  // namespace canoe
